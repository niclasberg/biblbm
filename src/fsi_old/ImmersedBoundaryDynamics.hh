#ifndef IMMERSEDBOUNDARYDYNAMICS_HH_
#define IMMERSEDBOUNDARYDYNAMICS_HH_
#include "ImmersedBoundaryDynamics.h"
#include "atomicBlock/dataProcessor3D.h"
#include "core/geometry3D.h"
#include <algorithm>
#include "CommunicationBuffer.hh"
#include "ParticleShapeFactory.h"
#include <sstream>
#include "MacroProcessors.h"
#include "Boundary.h"
#include <stdexcept>
#include "IO.h"
#include "Buffer.h"

#define FSI_PROFILE
#include "Profile.h"


namespace plb {

namespace fsi {

/*------ SolidNodeList ------*/
template<class T>
void SolidNodeList<T>::push_static_node(const SolidNode<T> & node)
{
	if(solid_node_capacity <= solid_node_count) {
		solid_node_capacity *= 2;
		SolidNode<T> * tmp = new SolidNode<T>[solid_node_capacity];
		std::copy(solid_nodes, solid_nodes+solid_node_count, tmp);
		delete [] solid_nodes;
		solid_nodes = tmp;
	}
	solid_nodes[solid_node_count++] = node;
}


/*------ ImmersedBoundaryDynamics3D ------*/
template<class T, template<typename U> class Descriptor, class Periodicity>
ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::ImmersedBoundaryDynamics3D(
		const MultiBlockLattice3D<T, Descriptor> & lattice,
		const ParticleShapeLibrary<T> & shape_library_)
: shape_library(shape_library_),
  management(lattice.getMultiBlockManagement()),
  num_particles(0),
  particle_max_radius(shape_library_.get_max_particle_radius()),
  arithmetic(Periodicity::create_arithmetic(lattice.getBoundingBox())),
  particle_force_communicator(particles, Periodicity::create_arithmetic(lattice.getBoundingBox())),
  particle_state_communicator(Periodicity::create_arithmetic(lattice.getBoundingBox())),
  fsi_force_communicator(Periodicity::create_arithmetic(lattice.getBoundingBox())),
  output_stream(0),
  fsi_subiterations(1),
  boundary(0),
  global_bounding_box(lattice.getBoundingBox()),
  nx(lattice.getNx()),
  ny(lattice.getNy()),
  nz(lattice.getNz()),
  is_init(false),
  acceleration_factor(0.45)
{
	initialize_domains();
	initialize_periodicity(lattice);
	initialize_fsi_communicator();
	set_interaction_cutoff_distance(1.0);
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::initialize_domains()
{
	// Get the block ids corresponding to the current MPI process
	std::vector<plint> local_blocks = management.getSparseBlockStructure().getLocalBlocks(management.getThreadAttribution());

	if(local_blocks.size() != 1) {
		pcerr << "The fsi module can only handle multiblock representations with 1 block per thread" << std::endl;
		exit(-1);
	}

	// Save local domain
	management.getSparseBlockStructure().getBulk(local_blocks[0], domain);

	// Save domain with envelope
	domain_with_fsi_envelope = domain;
	domain_with_fsi_envelope.enlarge_inplace(Dirac::half_support);

	// Bulk domain (nodes in this subdomain does not depend on entities in other processors)
	// By doing this subdivision, a lot of if-else statements can be avoided in the
	// cpu-intensive loops and branching can thus be kept at a minimum.
	bulk_domain = domain;
	bulk_domain.enlarge_inplace(-Dirac::half_support + 1.0);
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::initialize_periodicity(const MultiBlockLattice3D<T, Descriptor> & lattice)
{
	// Verify that the lattice has the same periodicity as that defined by the Periodicity template parameter
	if(lattice.periodicity().get(0) != Periodicity::get_x() ||
		lattice.periodicity().get(1) != Periodicity::get_y() ||
		lattice.periodicity().get(2) != Periodicity::get_z())
	{
		std::cerr << "The lattice must have the same periodicity as the fsi module." << std::endl
				<< "Check the 3rd template parameter of the ImmersedBoundaryDynamics3D class" << std::endl;
		exit(-1);
	}

	// Check if the local domain has edges that are connected by periodic boundary conditions
	has_periodic_edge[0] = has_periodic_edge[1] = has_periodic_edge[2] = false;
	if(Periodicity::get_x() && domain.x0 == global_bounding_box.x0 && domain.x1 == global_bounding_box.x1)
		has_periodic_edge[0] = true;
	if(Periodicity::get_y() && domain.y0 == global_bounding_box.y0 && domain.y1 == global_bounding_box.y1)
		has_periodic_edge[1] = true;
	if(Periodicity::get_z() && domain.z0 == global_bounding_box.z0 && domain.z1 == global_bounding_box.z1)
		has_periodic_edge[2] = true;
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::initialize_fsi_communicator()
{
	std::vector<plint> neighboring_processors;
	std::vector<std::pair<geo::Rect<T>, plint> > domain_proc_map;

	// Setup the communicator for the evaluation of the fsi forces
	// acting on the particle nodes
	create_proc_map(domain_with_fsi_envelope,
			Dirac::half_support,
			neighboring_processors,
			domain_proc_map);

	fsi_force_communicator.set_proc_list(neighboring_processors);
	fsi_force_communicator.set_domain_proc_map(domain_proc_map);
}

template<typename T, template< typename U > class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::create_proc_map(
		const geo::Rect<T> & search_domain,
		T envelope,
		std::vector<plint> & proc_list,
		std::vector<std::pair<geo::Rect<T>, plint> > & proc_domain_map) const
{
	proc_list.clear();
	proc_domain_map.clear();

	const SparseBlockStructure3D & block_structure = management.getSparseBlockStructure();
	const ThreadAttribution & thread_info = management.getThreadAttribution();

	// Find neighboring blocks
	std::vector<plint> neighboring_blocks;
	find_blocks_in_domain(search_domain, neighboring_blocks);

	// Get the processor ids corresponding to the neighboring blocks
	for(std::vector<plint>::iterator it = neighboring_blocks.begin(); it != neighboring_blocks.end(); ++it) {
		if(thread_info.isLocal(*it))
			continue;

		// Get mpi proc id
		plint proc_id = thread_info.getMpiProcess(*it);

		// Get the corresponding domain and add an envelope
		Box3D box;
		block_structure.getBulk(*it, box);
		geo::Rect<T> box_with_envelope(box);
		box_with_envelope.enlarge_inplace(envelope);

		proc_domain_map.push_back(std::make_pair(box_with_envelope, proc_id));

		// Avoid adding duplicates
		if(std::find(proc_list.begin(), proc_list.end(), proc_id) == proc_list.end())
			proc_list.push_back(proc_id);
	}
}

template<typename T, template< typename U > class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::find_blocks_in_domain(
		const geo::Rect<T> & bb,
		std::vector<plint> & neighboring_blocks
) const
{
	neighboring_blocks.clear();
	const std::map<plint,Box3D> & blocks = management.getSparseBlockStructure().getBulks();
	for(std::map<plint,Box3D>::const_iterator it = blocks.begin(); it != blocks.end(); ++it) {
		geo::Rect<T> domain_to_test(it->second);
		if(geo::does_intersect(bb, domain_to_test, arithmetic))
			neighboring_blocks.push_back(it->first);
	}
}

template<typename T, template< typename U > class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::set_interaction_cutoff_distance(T val)
{
	// Save interaction distance and its square
	interaction_distance = val;
	interaction_distance_sqr = val*val;

	// Update the interaction domain
	interaction_domain = domain;
	interaction_domain.enlarge_inplace(0.5 + interaction_distance);

	// Set the bounding domain
	domain_bounds = geo::Rect<T>(
		std::min(interaction_domain.x0, domain_with_fsi_envelope.x0),
		std::max(interaction_domain.x1, domain_with_fsi_envelope.x1),
		std::min(interaction_domain.y0, domain_with_fsi_envelope.y0),
		std::max(interaction_domain.y1, domain_with_fsi_envelope.y1),
		std::min(interaction_domain.z0, domain_with_fsi_envelope.z0),
		std::max(interaction_domain.z1, domain_with_fsi_envelope.z1)
	);

	rebuild_neighbor_proc_list();
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::rebuild_neighbor_proc_list()
{
	geo::Rect<T> test_domain = domain_bounds.enlarge(particle_max_radius);
	T envelope = std::max(Dirac::half_support, interaction_distance+0.5);

	std::vector<plint> neighboring_processors;
	std::vector<std::pair<geo::Rect<T>, plint> > domain_proc_map;

	create_proc_map(test_domain,
			envelope,
			neighboring_processors,
			domain_proc_map);

	particle_state_communicator.set_proc_list(neighboring_processors);
	particle_state_communicator.set_domain_proc_map(domain_proc_map);

	// Now setup the communicator for reduction of particle forces and torques
	create_proc_map(
			test_domain,
			envelope,
			neighboring_processors,
			domain_proc_map);

	particle_force_communicator.set_proc_list(neighboring_processors);
	particle_force_communicator.set_domain_proc_map(domain_proc_map);
}

/* Destructor */
template<class T, template<typename U> class Descriptor, class Periodicity>
ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::~ImmersedBoundaryDynamics3D()
{
	for(typename ObjMapType::iterator it = particles.begin(); it != particles.end(); ++it)
		delete it->second;
	if(output_stream)
		delete output_stream;
}

/* Dry iteration
 * Performs an iteration without taking fsi into account
 * */
template<class T, template<typename U> class Descriptor, class Periodicity>
	template<class InteractionFunctional>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::do_dry_iteration(
		const InteractionFunctional & interaction,
		T damping)
{
	if(this->iteration == 0) {
		// For the first iteration, no previous interaction forces has been computed and stored
		compute_interaction_forces(interaction);
	}

	// Update orientation and center of mass position of the particles
	move_particles(damping);

	// Compute interaction forces
	compute_interaction_forces(interaction);

	// Integrate momentum equations
	for(typename ObjMapType::iterator it = particles.begin(), it_end = particles.end(); it != it_end; ++it) {
		it->second->integrate_momentum_no_fsi(iteration, damping);
	}
}

/* Interaction force computation */
template<typename T, template< typename U > class Descriptor, class Periodicity>
	template<class InteractionFunctional>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::compute_interaction_forces(const InteractionFunctional & interaction)
{
	Profile::start_timer("Interaction force computation");

	plint num_inter = 0;

	// Reset collision forces and torques to zero and store the previous values
	for(typename ObjMapType::iterator it = particles.begin(); it != particles.end(); ++it) {
		it->second->coll_force_last = it->second->coll_force;
		it->second->coll_torque_last = it->second->coll_torque;
		it->second->coll_force.resetToZero();
		it->second->coll_torque.resetToZero();
	}

	// Particle-particle interactions
	if(particles.size() > 1) {
		for(typename ObjMapType::iterator it = particles.begin(); it != particles.end(); ++it) {
			// Only compute forces if the first particle is local
			// This is a small trick to ensure that the interaction force is not computed twice
			if(it->second->get_proc_id() == global::mpi().getRank()) {
				// Access the iterator to the next particle
				typename ObjMapType::iterator it2 = it;
				++it2;

				for(; it2 != particles.end(); ++it2) {
					bool ret = it->second->compute_collision_forces(*(it2->second), interaction, arithmetic);
					if(ret) ++num_inter;
				}
			}
		}
	}

	particle_force_communicator.pack_from_particles(ParticleForceCommunicatorType::CollForce);
	particle_force_communicator.send_and_receive_no_wait();

	//if(num_inter != 0)
	//	std::cout << global::mpi().getRank() << ": " << num_inter << " particle-particle interactions" << std::endl;

	// Particle-wall interactions
	if(has_boundary()) {
		for(typename ObjMapType::iterator it = particles.begin(); it != particles.end(); ++it) {
			it->second->compute_wall_collision_forces(*boundary, interaction);
		}
	}

	particle_force_communicator.finalize_send_and_receive();
	particle_force_communicator.unpack_to_particles(ParticleForceCommunicatorType::CollForce);

	Profile::stop_timer("Interaction force computation");
}

namespace detail {

template<class T>
void dump_particle_info(const RigidParticle3D<T> & p)
{
	// Print debug information
	std::cerr << "id: " << p.get_id() << std::endl
			  << "Position: (" << p.get_position()[0] << ", "
				  << p.get_position()[1] << ", "
				  << p.get_position()[2] << ")" << std::endl
			  << "Velocity: (" << p.get_velocity()[0] << ", "
				  << p.get_velocity()[1] << ", "
				  << p.get_velocity()[2] << ")" << std::endl
			  << "Orientation: (" << p.get_orientation()[0] << ", "
				  << p.get_orientation()[1] << ", "
				  << p.get_orientation()[2] << ", "
				  << p.get_orientation()[3] << ")" << std::endl
			  << "Angular velocity: (" << p.get_angular_velocity()[0] << ", "
				  << p.get_angular_velocity()[1] << ", "
				  << p.get_angular_velocity()[2] << ")" << std::endl
			  << "Fsi force: (" << p.get_force()[0] << ", "
				  << p.get_force()[1] << ", "
				  << p.get_force()[2] << ")" << std::endl
			  << "Fsi force (last timestep): (" << p.fsi_force_tmp[0] << ", "
				  << p.fsi_force_tmp[1] << ", "
				  << p.fsi_force_tmp[2] << ")" << std::endl
			  << "Fsi torque: (" << p.get_torque()[0] << ", "
				  << p.get_torque()[1] << ", "
				  << p.get_torque()[2] << ")" << std::endl
			  << "Fsi torque (last timestep): (" << p.fsi_torque_tmp[0] << ", "
				  << p.fsi_torque_tmp[1] << ", "
				  << p.fsi_torque_tmp[2] << ")" << std::endl
			  << "Collision force: (" << p.coll_force[0] << ", "
				  << p.coll_force[1] << ", "
				  << p.coll_force[2] << ")" << std::endl
			  << "Collision force (last timestep): (" << p.coll_force_tmp[0] << ", "
				  << p.coll_force_tmp[1] << ", "
				  << p.coll_force_tmp[2] << ")" << std::endl
			  << "Collision torque: (" << p.coll_torque[0] << ", "
				  << p.coll_torque[1] << ", "
				  << p.coll_torque[2] << ")" << std::endl
			  << "Collision torque (last timestep): (" << p.coll_torque_tmp[0] << ", "
				  << p.coll_torque_tmp[1] << ", "
				  << p.coll_torque_tmp[2] << ")" << std::endl;
}

}

/* Particle movement */
template<typename T, template< typename U > class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::move_particles(T damping)
{
	Profile::start_timer("Move particles");

	ObjMapType updated_particle_map;
	particle_state_communicator.clear_send_buffer();

	// Move all local particles
	for(typename ObjMapType::iterator it = particles.begin(), it_end = particles.end(); it != it_end; ++it) {
		if(it->second->get_proc_id() == global::mpi().getRank()) {
			it->second->move(damping);

			// Map position to periodic grid
			arithmetic.remap_position(it->second->get_position());

			// Pack for neighboring processors
			particle_state_communicator.pack_particle(*it->second);

			// Insert into the updated map if the particle still intersects this processor's domain
			if(geo::does_intersect(it->second->get_bounding_sphere(), domain_bounds, arithmetic))
				updated_particle_map[it->first] = it->second;
			else {
				// Out of bounds check
				if(management.getSparseBlockStructure().
						locate(std::floor(it->second->get_position()[0]+0.5), std::floor(it->second->get_position()[1]+0.5), std::floor(it->second->get_position()[2]+0.5))
						== -1) {
					std::cerr << "Warning: A particle is out of bounds!" << std::endl;
					detail::dump_particle_info(*(it->second));
					throw std::runtime_error("A particle is out of range");
				}
			}
		}
	}

	// Send and receive the new state of the particles (non-blocking)
	particle_state_communicator.send_and_receive_no_wait();

	// Update the solid nodes of the local particles
	local_node_list.clear();
	envelope_node_list.clear();

	for(typename ObjMapType::iterator it = updated_particle_map.begin(),
			it_end = updated_particle_map.end(); it != it_end; ++it) {
		add_solid_nodes_from_particle(*it->second);

		// Remove the nodes from the old object_map
		particles.erase(it->first);
	}

	// Receive the new state of the non-local particles
	particle_state_communicator.finalize_send_and_receive();

	// Update the solid nodes of the non-local particles
	for(typename ParticleStateCommunicatorType::DataIterator it = particle_state_communicator.data_begin(),
			it_max = particle_state_communicator.data_end(); it != it_max; ++it) {

		// Check if the particle already exists in the domain
		ObjMapIterator p_it = particles.find(it->obj_id);
		RigidParticle3D<T> * particle;
		if(p_it != particles.end()) {
			particle = p_it->second;
			particles.erase(p_it);
		} else {
			particle = new RigidParticle3D<T>(shape_library.get_by_id(it->shape_id));
			particle->set_id(it->obj_id);
		}

		// Copy received data
		particle->get_position() = it->pos;
		particle->get_angular_velocity() = it->ang_vel;
		particle->get_velocity() = it->vel;
		particle->get_orientation() = it->orientation;
		particle->get_density() = it->density;
		particle->get_scale() = it->scale;
		particle->get_force() = it->force;
		particle->get_torque() = it->torque;
		particle->coll_force = it->force_coll;
		particle->coll_force_last = it->force_coll_last;
		particle->coll_torque = it->torque_coll;
		particle->coll_torque_last = it->torque_coll_last;
		particle->ang_momentum = it->ang_momentum;
		particle->update_rotation_matrix();

		updated_particle_map[it->obj_id] = particle;
		add_solid_nodes_from_particle(*particle);
	}

	// Remove all particles no longer residing in the domain
	for(typename ObjMapType::iterator it = particles.begin(); it != particles.end(); ++it)
		delete it->second;

	particles.swap(updated_particle_map);

	// Set particle processor ids
	for(typename ObjMapType::iterator it = particles.begin(); it != particles.end(); ++it)
		set_particle_proc_id(*it->second);

	Profile::stop_timer("Move particles");
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::set_particle_proc_id(RigidParticle3D<T> & p)
{
	plint particle_block_id = management.getSparseBlockStructure().
				locate(std::floor(p.get_position()[0]+0.5), std::floor(p.get_position()[1]+0.5), std::floor(p.get_position()[2]+0.5));

	if(particle_block_id == -1) {
		std::cerr << "Warning: A particle is out of bounds!" << std::endl;
		detail::dump_particle_info(p);
		throw std::runtime_error("A particle is out of range");
	}

	p.set_proc_id(management.getThreadAttribution().getMpiProcess(particle_block_id));
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::add_solid_nodes_from_particle(RigidParticle3D<T> & p)
{
	// First do an overlap test to eliminate particles that are not overlapping this domain
	if( ! geo::does_intersect(p.get_bounding_sphere(), domain_with_fsi_envelope, arithmetic))
		return;

	Array<T, 3> pos;
	SolidNode<T> node;
	ParticleNodeInfo<T> node_info;
	const Array<T, 3> mid(0.5 * (domain.x0 + domain.x1),
						  0.5 * (domain.y0 + domain.y1),
						  0.5 * (domain.z0 + domain.z1));

	// Generate nodes
	for(pluint i = 0, i_max = p.get_vertex_count(); i < i_max; ++i) {
		p.get_centroid_world_frame(i, pos);

		// If the grid is periodic, move the position a multiple of the global domain extent so that
		// the position is as close to the center of the local domain as possible
		arithmetic.shift_periodically_to_minimize_distance_to(mid, pos);

		// Check if the node is contained in this domain
		if(domain_with_fsi_envelope.contains_or_on_boundary(pos))
		{
			// Generate node
			p.get_node_info(i, node_info);
			node.pos = pos;
			node.vel = node_info.vel;
			node.area = node_info.area;
			node.pos_rel = node_info.pos_rel;
			node.normal = node_info.normal;
			node.obj_id = p.get_id();
			node.node_id = i;
			node.particle = &p;

			// Check if the point is in the bulk or envelope domain
			if(bulk_domain.contains(pos)) {
				local_node_list.push_node(node);
			} else {
				envelope_node_list.push_node(node);
			}
		}
	}
	//std::cout << "Proc " << global::mpi().getRank() << ", envelope node count = " << envelope_node_list.size() << ", local node count = " << local_node_list.size() << std::endl;
}

/* Fsi iteration */
template<class T, template<typename U> class Descriptor, class Periodicity>
	template<class InteractionFunctional>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::do_fsi(
		Box3D dom,
		BlockLattice3D<T, Descriptor> & lattice,
		ScalarField3D<T> & density,
		TensorField3D<T, 3> & velocity,
		const Array<T, 3> & external_force,
		const InteractionFunctional & interaction)
{
	if(is_init) {
		// Update orientation and center of mass position of the particles
		move_particles();

		// Compute interaction forces for the new configuration
		compute_interaction_forces(interaction);

		// Predict the new linear and angular velocities
		for(typename ObjMapType::iterator it = particles.begin(), it_end = particles.end(); it != it_end; ++it) {
			it->second->integrate_momentum_1st_stage(iteration, external_force);
			it->second->get_force().resetToZero();
			it->second->get_torque().resetToZero();
			//std::cout << "0->" << it->second->get_angular_velocity()[2] << std::endl;
		}
		update_solid_nodes_velocities();

		for(pluint iter = 0; iter < fsi_subiterations; ++iter) {
			compute_fsi_forces(dom, lattice, density, velocity);

			for(typename ObjMapType::iterator it = particles.begin(), it_end = particles.end(); it != it_end; ++it) {
				it->second->integrate_momentum_2nd_stage(iter);
				it->second->get_torque() += it->second->fsi_torque_tmp;
				it->second->get_force() += it->second->fsi_force_tmp;
			}
			update_solid_nodes_velocities();
		}

	} else {
		compute_interaction_forces(interaction);

		for(pluint it = 0; it < fsi_subiterations; ++it) {
			compute_fsi_forces(dom, lattice, density, velocity);
			for(typename ObjMapType::iterator it = particles.begin(), it_end = particles.end(); it != it_end; ++it) {
				it->second->get_torque() += it->second->fsi_torque_tmp;
				it->second->get_force() += it->second->fsi_force_tmp;
			}
		}

		is_init = true;
	}
}

inline void restrict_box(Box3D & box, const Box3D & domain)
{
	if(box.x0 < domain.x0) box.x0 = domain.x0;
	if(box.x1 > domain.x1) box.x1 = domain.x1;
	if(box.y0 < domain.y0) box.y0 = domain.y0;
	if(box.y1 > domain.y1) box.y1 = domain.y1;
	if(box.z0 < domain.z0) box.z0 = domain.z0;
	if(box.z1 > domain.z1) box.z1 = domain.z1;
}

#ifdef FSI_COMPUTE_RESIDUAL
	double delta_u_max = 0;
	double delta_u_l2 = 0;
#endif

/*
 * Fsi force computation
 */
template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::compute_fsi_forces(
		Box3D dom,
		BlockLattice3D<T, Descriptor> & lattice,
		ScalarField3D<T> & density,
		TensorField3D<T, 3> & velocity)
{
	// Reset temporary forces on particles
	for(typename ObjMapType::iterator it = particles.begin(), it_end = particles.end(); it != it_end; ++it) {
		it->second->fsi_force_tmp.resetToZero();
		it->second->fsi_torque_tmp.resetToZero();
	}

	// Process nodes in the fsi envelope
	compute_envelope_fsi_forces(dom, lattice, density, velocity);

	// Send and receive (non-blocking)
	fsi_force_communicator.send_and_receive_no_wait();

	// Process bulk nodes while waiting for the data from the neighboring processors
	compute_bulk_fsi_forces(dom, lattice, density, velocity);

	// Synchronize with neighboring processors
	Profile::start_timer("Fsi force synchronization");
	fsi_force_communicator.finalize_send_and_receive();
	Profile::stop_timer("Fsi force synchronization");

	// Spread the fsi force contribution from neighboring processors
	spread_envelope_fsi_forces(dom, lattice, density, velocity);

	// Communicate particle forces
	particle_force_communicator.pack_from_particles(ParticleForceCommunicatorType::TmpForce);
	particle_force_communicator.send_and_receive_no_wait();

	// Update the fluid velocity in the domain
	update_fluid_velocity(dom, lattice, velocity);

	// Wait for data to arrive
	Profile::start_timer("Particle force synchronization");
	particle_force_communicator.finalize_send_and_receive();
	particle_force_communicator.unpack_to_particles(ParticleForceCommunicatorType::TmpForce);
	Profile::stop_timer("Particle force synchronization");

#ifdef FSI_COMPUTE_RESIDUAL
	T delta_u_max_recv, delta_u_l2_recv;
	global::mpi().reduce(delta_u_max, delta_u_max_recv, MPI_MAX);
	global::mpi().reduce(delta_u_l2, delta_u_l2_recv, MPI_SUM);

	if(global::mpi().isMainProcessor())
		std::cout << std::scientific << delta_u_max_recv/acceleration_factor << " " << std::sqrt(delta_u_l2_recv)/acceleration_factor << std::endl;
	std::cout.unsetf(std::ios_base::floatfield);
#endif
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::compute_envelope_fsi_forces(
		Box3D dom,
		BlockLattice3D<T, Descriptor> & lattice,
		ScalarField3D<T> & density,
		TensorField3D<T, 3> & velocity)
{
	Profile::start_timer("Envelope interpolation");

	// Offsets between the lattices
	Dot3D density_offset = computeRelativeDisplacement(lattice, density);
	Dot3D velocity_offset = computeRelativeDisplacement(lattice, velocity);
	Dot3D offset = lattice.getLocation();

	// Temporary variables
	Array<T, 3> vel_fluid;
	Box3D d_bounds;
	T dirac_x, dirac_xy, dirac_val;
	std::vector<T> dirac_values(Dirac::support * Dirac::support * Dirac::support);

	fsi_force_communicator.clear_data();
	for(pluint ni = 0, ni_max = envelope_node_list.size(); ni < ni_max; ++ni) {
		SolidNode<T> & node = envelope_node_list[ni];

		const Array<T, 3> & pos = node.pos;
		const Array<T, 3> pos_rel(pos[0] - offset.x, pos[1] - offset.y, pos[2] - offset.z);

		// Find compact support region of the Dirac function
		get_dirac_compact_support_box<T, Dirac>(pos_rel, d_bounds);

		// Interpolate fluid velocity
		vel_fluid.resetToZero();
		pluint cnt = 0;
		for(plint i = d_bounds.x0; i <= d_bounds.x1; ++i) {
			plint i2 = i;
			if(has_periodic_edge[0]) arithmetic.remap_index_x(i2, offset.x);
			if(i2 < dom.x0 || i2 > dom.x1) continue;

			dirac_x = sampled_dirac.eval((T)i - pos_rel[0]);
			for(plint j = d_bounds.y0; j <= d_bounds.y1; ++j) {
				plint j2 = j;
				if(has_periodic_edge[1]) arithmetic.remap_index_y(j2, offset.y);
				if(j2 < dom.y0 || j2 > dom.y1) continue;

				dirac_xy = dirac_x * sampled_dirac.eval((T)j - pos_rel[1]);
				for(plint k = d_bounds.z0; k <= d_bounds.z1; ++k) {
					plint k2 = k;
					if(has_periodic_edge[2]) arithmetic.remap_index_z(k2, offset.z);
					if(k2 < dom.z0 || k2 > dom.z1) continue;

					dirac_val = dirac_xy * sampled_dirac.eval((T)k - pos_rel[2]);
					vel_fluid += (dirac_val * velocity.get(i2 + velocity_offset.x,
												 j2 + velocity_offset.y,
												 k2 + velocity_offset.z));
					dirac_values[cnt] = dirac_val;
					++cnt;
				}
			}
		}

		// The node velocity is weighted by #nodes in processor / #nodes in dirac support region. This is a
		// hack so that the total force acting on the fluid adds up to the correct value from the contributions
		// from different processors.
		const Array<T, 3> delta_u = ((T) cnt / (T)(Dirac::support*Dirac::support*Dirac::support) * node.vel - vel_fluid)
												* node.area * acceleration_factor;

		// Spread velocity correction to the fluid
		T factor = 0;
		cnt = 0;
		for(plint i = d_bounds.x0; i <= d_bounds.x1; ++i) {
			plint i2 = i;
			if(has_periodic_edge[0]) arithmetic.remap_index_x(i2, offset.x);
			if(i2 < dom.x0 || i2 > dom.x1) continue;
			for(plint j = d_bounds.y0; j <= d_bounds.y1; ++j) {
				plint j2 = j;
				if(has_periodic_edge[1]) arithmetic.remap_index_y(j2, offset.y);
				if(j2 < dom.y0 || j2 > dom.y1) continue;
				for(plint k = d_bounds.z0; k <= d_bounds.z1; ++k) {
					plint k2 = k;
					if(has_periodic_edge[2]) arithmetic.remap_index_z(k2, offset.z);
					if(k2 < dom.z0 || k2 > dom.z1) continue;

					const T multiple = density.get(i2 + density_offset.x,
												 j2 + density_offset.y,
												 k2 + density_offset.z)
											 * dirac_values[cnt] * 2.0;
					factor += multiple;

					T * g = lattice.get(i2, j2, k2).getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
					g[0] += delta_u[0] * multiple;
					g[1] += delta_u[1] * multiple;
					g[2] += delta_u[2] * multiple;

					++cnt;
				}
			}
		}

		// Add force contribution to the appropriate particle
		const Array<T, 3> force_local = -factor * delta_u;
		node.particle->fsi_force_tmp += force_local;
		node.particle->fsi_torque_tmp += crossProduct(node.pos_rel, force_local);

		// Pack for mpi send
		const VelocityInterpolationNode<T> send_node = {delta_u, node.pos, node.pos_rel, node.obj_id};
		fsi_force_communicator.add_data(send_node);
	}

	Profile::stop_timer("Envelope interpolation");
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::compute_bulk_fsi_forces(
		Box3D dom,
		BlockLattice3D<T, Descriptor> & lattice,
		ScalarField3D<T> & density,
		TensorField3D<T, 3> & velocity)
{
	Profile::start_timer("Bulk fsi");

	// Offsets between the lattices
	Dot3D density_offset = computeRelativeDisplacement(lattice, density);
	Dot3D velocity_offset = computeRelativeDisplacement(lattice, velocity);
	Dot3D offset = lattice.getLocation();

	// Temporary variables
	Array<T, 3> vel_fluid;
	Box3D d_bounds;
	T dirac_x, dirac_xy, dirac_val;
	std::vector<T> dirac_values(Dirac::support * Dirac::support * Dirac::support);

	for(pluint ni = 0, ni_max = local_node_list.size(); ni < ni_max; ++ni) {
		SolidNode<T> & node = local_node_list[ni];
		const Array<T, 3> & pos = node.pos;
		const Array<T, 3> pos_rel(pos[0] - offset.x, pos[1] - offset.y, pos[2] - offset.z);

		// Find compact support region of the Dirac function
		get_dirac_compact_support_box<T, Dirac>(pos_rel, d_bounds);

		// Interpolate fluid velocity
		// Since these nodes are completely local, no special treatment such as
		// out of bounds checks or periodicity imposing is required.
		vel_fluid.resetToZero();
		plint cnt = 0;
		for(plint i = d_bounds.x0; i <= d_bounds.x1; ++i) {
			dirac_x = sampled_dirac.eval((T)i - pos_rel[0]);
			for(plint j = d_bounds.y0; j <= d_bounds.y1; ++j) {
				dirac_xy = dirac_x * sampled_dirac.eval((T)j - pos_rel[1]);
				for(plint k = d_bounds.z0; k <= d_bounds.z1; ++k) {
					dirac_val = dirac_xy * sampled_dirac.eval((T)k - pos_rel[2]);
					dirac_values[cnt] = dirac_val;
					vel_fluid += dirac_val * velocity.get(i + velocity_offset.x,
													 j + velocity_offset.y,
													 k + velocity_offset.z);
					++cnt;
				}
			}
		}

		// Area-weighted velocity correction
		const Array<T, 3> delta_u = (node.vel - vel_fluid) * node.area * acceleration_factor;

#ifdef FSI_COMPUTE_RESIDUAL
		for(pluint k = 0; k < 3; ++k)
			delta_u_max = std::abs(delta_u[k]) > delta_u_max ? std::abs(delta_u[k]) : delta_u_max;
		delta_u_l2 += normSqr(delta_u);
#endif

		// Spread the velocity correction and compute the corresponding force addition
		cnt = 0;
		T factor = 0;
		for(plint i = d_bounds.x0; i <= d_bounds.x1; ++i)
			for(plint j = d_bounds.y0; j <= d_bounds.y1; ++j)
				for(plint k = d_bounds.z0; k <= d_bounds.z1; ++k) {
					const T multiple = density.get(i + density_offset.x,
										 j + density_offset.y,
										 k + density_offset.z)
												 * dirac_values[cnt] * 2.0;
					factor += multiple;
					T * g = lattice.get(i, j, k).getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
					g[0] += delta_u[0] * multiple;
					g[1] += delta_u[1] * multiple;
					g[2] += delta_u[2] * multiple;

					++cnt;
				}

		// Add force and torque contribution to particle
		node.particle->fsi_force_tmp -= factor * delta_u;
		node.particle->fsi_torque_tmp -= factor * crossProduct(node.pos_rel, delta_u);
	}
	Profile::stop_timer("Bulk fsi");
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::spread_envelope_fsi_forces(
		Box3D dom,
		BlockLattice3D<T, Descriptor> & lattice,
		ScalarField3D<T> & density,
		TensorField3D<T, 3> & velocity)
{
	Profile::start_timer("Fsi force envelope spreading");

	// Offsets between the lattices
	Dot3D density_offset = computeRelativeDisplacement(lattice, density);
	Dot3D velocity_offset = computeRelativeDisplacement(lattice, velocity);
	Dot3D offset = lattice.getLocation();

	// Temporary variables
	Array<T, 3> vel_fluid;
	Box3D d_bounds;
	T dirac_x, dirac_xy, dirac_val;

	const Array<T, 3> mid(0.5 * (domain.x0 + domain.x1),
							  0.5 * (domain.y0 + domain.y1),
							  0.5 * (domain.z0 + domain.z1));
	for(typename FsiForceCommunicatorType::IteratorType it = fsi_force_communicator.data_begin(),
			it_end = fsi_force_communicator.data_end();
			it != it_end; ++it) {

		Array<T, 3> pos = it->pos;
		arithmetic.shift_periodically_to_minimize_distance_to(mid, pos);
		const Array<T, 3> pos_rel(
				pos[0] - offset.x,
				pos[1] - offset.y,
				pos[2] - offset.z);

		get_dirac_compact_support_box<T, Dirac>(pos_rel, d_bounds);

		// Spread velocity correction to the fluid
		T factor = 0;
		for(plint i = d_bounds.x0; i <= d_bounds.x1; ++i) {
			plint i2 = i;
			if(has_periodic_edge[0]) arithmetic.remap_index_x(i2, offset.x);
			if(i2 < dom.x0 || i2 > dom.x1) continue;

			dirac_x = sampled_dirac.eval((T)i - pos_rel[0]);
			for(plint j = d_bounds.y0; j <= d_bounds.y1; ++j) {
				plint j2 = j;
				if(has_periodic_edge[1]) arithmetic.remap_index_y(j2, offset.y);
				if(j2 < dom.y0 || j2 > dom.y1) continue;

				dirac_xy = dirac_x * sampled_dirac.eval((T)j - pos_rel[1]);
				for(plint k = d_bounds.z0; k <= d_bounds.z1; ++k) {
					plint k2 = k;
					if(has_periodic_edge[2]) arithmetic.remap_index_z(k2, offset.z);
					if(k2 < dom.z0 || k2 > dom.z1) continue;

					const T multiple = density.get(i2 + density_offset.x,
												 j2 + density_offset.y,
												 k2 + density_offset.z)
											 * dirac_xy * sampled_dirac.eval((T)k - pos_rel[2]) * 2.0;
					factor += multiple;

					T * g = lattice.get(i2, j2, k2).getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
					g[0] += it->delta_u[0] * multiple;
					g[1] += it->delta_u[1] * multiple;
					g[2] += it->delta_u[2] * multiple;
				}
			}
		}

		// Add force contribution to the appropriate particle
		RigidParticle3D<T> * particle = particles.at(it->obj_id);
		const Array<T, 3> force_local = -factor * it->delta_u;
		particle->fsi_force_tmp += force_local;
		particle->fsi_torque_tmp += crossProduct(it->pos_rel, force_local);
	}

	Profile::stop_timer("Fsi force envelope spreading");
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::update_fluid_velocity(
		Box3D dom,
		BlockLattice3D<T, Descriptor> & lattice,
		TensorField3D<T, 3> & velocity)
{
	Profile::start_timer("Update velocity");

	Dot3D velocity_offset = computeRelativeDisplacement(lattice, velocity);

	if(particles.size() != 0) {
		for(plint i = dom.x0; i <= dom.x1; ++i) {
			for(plint j = dom.y0; j <= dom.y1; ++j) {
				for(plint k = dom.z0; k <= dom.z1; ++k) {
					lattice.get(i, j, k).computeVelocity(
							velocity.get(i+velocity_offset.x,
									j+velocity_offset.y,
									k+velocity_offset.z));
				}
			}
		}
	}

	Profile::stop_timer("Update velocity");
}

template<typename T, template< typename U > class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::update_solid_nodes_velocities()
{
	// Update local nodes
	for(pluint i = 0, i_max = local_node_list.size(); i < i_max; ++i) {
		SolidNode<T> & node = local_node_list[i];
		node.particle->get_node_velocity_world_frame(node.node_id, node.vel);
	}

	// Update envelope nodes
	for(pluint i = 0, i_max = envelope_node_list.size(); i < i_max; ++i) {
		SolidNode<T> & node = envelope_node_list[i];
		node.particle->get_node_velocity_world_frame(node.node_id, node.vel);
	}
}

/*
 * add_particle
 * Adds a particle to the simulation. Particles not in the domain are discarded.
 */
template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::add_particle(RigidParticle3D<T> * particle)
{
	// Does the neighboring processor list need to be rebuilt?
	if(particle->get_radius() > particle_max_radius) {
		particle_max_radius = particle->get_radius();
		rebuild_neighbor_proc_list();
	}

	// Create an id for the particle.
	plint p_id = get_next_id();

	// Check if the particle's bounding sphere is in this domain
	if( ! geo::does_intersect(particle->get_bounding_sphere(), domain_bounds, arithmetic))
		return;

	// Copy the particle
	RigidParticle3D<T> * p_copy = particle->clone();

	// Set particle id
	p_copy->set_id(p_id);

	//std::cout << global::mpi().getRank() << "-> " << p_copy->get_id() << ": " << p_copy->get_position()[0] << ", " << p_copy->get_position()[1] << ", " << p_copy->get_position()[2] << std::endl;

	set_particle_proc_id(*p_copy);

	// Save the particle
	particles[p_id] = p_copy;
	add_solid_nodes_from_particle(*p_copy);
}

/*
 * Particle overlap testing
 */
template<class T, template<typename U> class Descriptor, class Periodicity>
bool ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::has_overlapping_particles() const
{
	int overlap_found = 0;
/*	for(typename ObjMapType::iterator it = particles.begin(); it != particles.end(); ++it) {
		for(typename ObjMapType::iterator it2 = it+1; it2 != particles.end(); ++it2) {
			// Do bounding sphere test
			if(does_intersect(it->second->get_bounding_sphere(), it2->second->get_bounding_sphere(), arithmetic)) {
				if(it->second->intersects(*(it2->second))) {
					overlap_found = 1;
					break;
				}
			}
		}

		if(overlap_found)
			break;
	}

	global::mpi().reduceAndBcast(overlap_found, MPI_LAND);
*/
	return overlap_found;
}

/*
 * Count total number of particles
 * NOTE: this method performs a reduction over all processors and is therefore SLOW. Use with care.
 */
template<class T, template<typename U> class Descriptor, class Periodicity>
pluint ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::count_particles() const
{
	plint num_particles = 0;
	for(typename ObjMapType::const_iterator it = particles.begin(); it != particles.end(); ++it) {
		if(it->second->get_proc_id() == global::mpi().getRank())
			++num_particles;
	}

#ifdef PLB_MPI_PARALLEL
	global::mpi().reduceAndBcast(num_particles, MPI_SUM);
#endif

	return num_particles;
}


/*
 * Output
 */
template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::write_particles_as_vtk(pluint iter)
{
	unsigned int p_id_width = std::floor(std::log10(num_particles-1)) + 1;

	for(typename ObjMapType::iterator it = particles.begin(); it != particles.end(); ++it) {
		if(it->second->get_proc_id() == global::mpi().getRank()) {
			std::stringstream fNameStream;
			fNameStream << global::directories().getVtkOutDir() << "obj" << std::setfill('0') << std::setw(p_id_width) << it->first
					<< "-" << std::setfill('0') << std::setw(7) << iter
					<<".vtu";

			std::ofstream out;
			out.open( fNameStream.str().c_str(), std::ios_base::trunc);

			if( ! out.good()) {
				std::cerr << "Could not open the file " << fNameStream.str() << " for writing." << std::endl;
			} else {
				it->second->write_to_stream_as_vtk(out);
				out.close();
			}
		}
	}
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::write_particle_states(pluint iter)
{
	if( ! output_stream) {
		// Create file name
		std::stringstream fNameStream;
		fNameStream << global::directories().getVtkOutDir()
							  << "objsP" << std::setfill('0') << std::setw(6) << global::mpi().getRank() <<".txt";
		output_stream = new std::ofstream(fNameStream.str().c_str(), std::ios_base::out);

		// Print header
		*output_stream << "iteration particle_id position_x position_y position_z velocity_x velocity_y velocity_z "
				<< "orientation_scalar orientation_x orientation_y orientation_z "
				<< "angular_velocity_x angular_velocity_y angular_velocity_z "
				<< "fsi_force_x fsi_force_y fsi_force_z fsi_torque_x fsi_torque_y fsi_torque_z "
				<< "coll_force_x coll_force_y coll_force_z coll_torque_x coll_torque_y coll_torque_z "
				<< "kinetic_energy" << std::endl;
	}

	for(typename ObjMapType::iterator it = particles.begin(); it != particles.end(); ++it) {
		if(it->second->get_proc_id() == global::mpi().getRank()) {
			RigidParticle3D<T> * p = it->second;
			*output_stream << iter << " "
							<< p->get_id() << " "
							<< p->get_position()[0] << " "
							<< p->get_position()[1] << " "
							<< p->get_position()[2] << " "
							<< p->get_velocity()[0] << " "
							<< p->get_velocity()[1] << " "
							<< p->get_velocity()[2] << " "
							<< p->get_orientation()[0] << " "
							<< p->get_orientation()[1] << " "
							<< p->get_orientation()[2] << " "
							<< p->get_orientation()[3] << " "
							<< p->get_angular_velocity()[0] << " "
							<< p->get_angular_velocity()[1] << " "
							<< p->get_angular_velocity()[2] << " "
							<< p->get_force()[0] << " "
							<< p->get_force()[1] << " "
							<< p->get_force()[2] << " "
							<< p->get_torque()[0] << " "
							<< p->get_torque()[1] << " "
							<< p->get_torque()[2] << " "
							<< p->coll_force[0] << " "
							<< p->coll_force[1] << " "
							<< p->coll_force[2] << " "
							<< p->coll_torque[0] << " "
							<< p->coll_torque[1] << " "
							<< p->coll_torque[2] << " "
							<< p->compute_kinetic_energy()
							<< std::endl;
		}
	}
}


template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::save_checkpoint(FileName file_name) const
{
	file_name.defaultPath(global::directories().getOutputDir());
	file_name.defaultExt("plb");

	// Output size of each particle
	const pluint particle_size = 10 * sizeof(Array<T, 3>) + sizeof(Quaternion<T>) + sizeof(pluint) + 2*sizeof(plint);
	Buffer<char> buffer(particle_size);

	// Create name for output folder
	std::string folder(file_name.getPath() + std::string("/") + file_name.getName());

	// Create folder (this mkdir implementation is thread safe)
	io::mkdir(folder.c_str());

	// Write header file
	if(global::mpi().isMainProcessor()) {
		// Open header file (the header file has the same name as the folder, but with extension .plb)
		std::ofstream out(file_name.get().c_str(), std::ios::out);

		// Write shared particle information
		out << particle_max_radius << std::endl;
		out << num_particles << std::endl;

		// Write file information
		out << global::mpi().getSize() << std::endl;
		for(plint i = 0; i < global::mpi().getSize(); ++i) {
			std::stringstream fss;
			fss << folder << "/" << std::setfill('0') << std::setw(6) << i << ".plb";
			out << fss.str() << std::endl;
		}
		out.close();
	}

	// Open output file
	std::stringstream ss;
	ss << folder << "/" << std::setfill('0') << std::setw(6) << global::mpi().getRank() << ".plb";
	std::ofstream out(ss.str().c_str(), std::ios::binary|std::ios::out);
	if( ! out.is_open()) {
		std::cerr << "Could not open the checkpoint file " << ss.str() << " for writing" << std::endl;
		return;
	}

	// Count local particles and write to file
	plint p_size = 0;
	for(typename ObjMapType::const_iterator it = particles.begin(); it != particles.end(); ++it)
		if(it->second->get_proc_id() == global::mpi().getRank())
			++p_size;

	out.write((char *) &p_size, sizeof(plint));

	// Serialize the particles
	for(typename ObjMapType::const_iterator it = particles.begin(); it != particles.end(); ++it) {
		if(it->second->get_proc_id() == global::mpi().getRank()) {
			buffer.rewind_ptr();
			buffer.pack(it->second->get_shape_id());
			buffer.pack(it->second->get_id());
			buffer.pack(it->second->get_position());
			buffer.pack(it->second->get_velocity());
			buffer.pack(it->second->get_angular_velocity());
			buffer.pack(it->second->get_orientation());
			buffer.pack(it->second->get_angular_velocity());
			buffer.pack(it->second->get_density());
			buffer.pack(it->second->get_force());
			buffer.pack(it->second->get_torque());
			buffer.pack(it->second->coll_force);
			buffer.pack(it->second->coll_force_last);
			buffer.pack(it->second->coll_torque);
			buffer.pack(it->second->coll_torque_last);

			buffer.rewind_ptr();
			out.write(buffer.get_data_ptr(), particle_size);
		}
	}
	out.close();
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::load_checkpoint(FileName file_name)
{
	file_name.defaultPath(global::directories().getInputDir());
	file_name.defaultExt("plb");

	// Clear everything particle-related
	particles.clear();
	local_node_list.clear();
	envelope_node_list.clear();

	// Communication
	const pluint particle_size = 10 * sizeof(Array<T, 3>) + sizeof(Quaternion<T>) + sizeof(pluint) + 2*sizeof(plint);
	const pluint particles_per_message = 40;
	const pluint message_envelope_size = sizeof(bool) + sizeof(plint);
	const pluint message_size = message_envelope_size + particles_per_message * particle_size;
	Buffer<char> buffer(message_size);
	pluint particles_received;

	// Temporal variables
	plint num_local_particles;
	plint p_id;
	pluint shape_id;

	std::vector<std::string> file_list;
	file_list.clear();
	T particle_max_radius_;
	plint num_particles_;

	// Read header file
	if(global::mpi().isMainProcessor()) {
		std::ifstream in(file_name.get().c_str(), std::ios::in);

		if( ! in.is_open()) {
			pcerr << "Could not open the checkpoint file " << file_name.get().c_str() << ". Does it exist?" << std::endl;
			exit(-1);
		}

		// Read number of particles and particle max radius
		in >> particle_max_radius_;
		in >> num_particles_;

		// Read file names
		pluint num_files;
		std::string fname;
		in >> num_files;
		while(std::getline(in, fname)) {
			if( ! fname.empty())
				file_list.push_back(fname);
		}
		in.close();
	}

	if(file_list.empty())
		pcerr << "Possible error in the checkpoint header file " << file_name.get() << ": No checkpoint files specified" << std::endl;

	global::mpi().bCast(&particle_max_radius_, 1);
	global::mpi().bCast(&num_particles_, 1);

	this->particle_max_radius = particle_max_radius_;
	this->num_particles = num_particles_;

	// Now that we know the maximal radius of all particles, rebuild the neighbor processor list
	rebuild_neighbor_proc_list();

	// Read all particle states. Note that the main processor does all the reading,
	// and the read information is then broadcasted to all other processors.
	bool keep_reading = true;
	std::vector<std::string>::iterator file_it = file_list.begin();
	plint particles_left_in_file = 0;
	std::ifstream * stream = 0;

	while(keep_reading) {
		if(global::mpi().isMainProcessor()) {
			// Pack dummy values for keep_reading and number of particles
			buffer.rewind_ptr();
			buffer.pack(keep_reading);
			buffer.pack(particles_received);

			// Read particles_per_message particle states from checkpoint files
			particles_received = 0;
			for(pluint i = 0; i < particles_per_message; ++i) {
				// If we have reached the end of the current file, proceed with the next one
				while(particles_left_in_file == 0) {
					if(stream)
						delete stream;

					if(file_it == file_list.end()) {
						keep_reading = false;
						break;
					}

					stream = new std::ifstream(file_it->c_str(), std::ios::binary | std::ios::in);

					if( ! stream->is_open())
						pcerr << "Could not open the checkpoint file " << *file_it << std::endl;

					// Read number of particles in file
					stream->read((char *) &particles_left_in_file, sizeof(plint));

					++file_it;
				}

				if( ! keep_reading)
					break;

				// Read particle data
				stream->read(buffer.get_data_ptr(), particle_size);
				buffer.advance_ptr(particle_size);

				++particles_received;
				--particles_left_in_file;
			}

			// Set appropriate values for keep_reading and particles_received
			buffer.rewind_ptr();
			buffer.pack(keep_reading);
			buffer.pack(particles_received);
		}

		// Broadcast the read data
		buffer.rewind_ptr();
		global::mpi().bCast(buffer.get_data_ptr(), message_size);

		buffer.unpack(keep_reading);
		buffer.unpack(particles_received);

		for(plint i = 0; i < particles_received; ++i) {
			// Read shape id and create the particle
			buffer.unpack(shape_id);
			RigidParticle3D<T> * particle = new RigidParticle3D<T>(shape_library.get_by_id(shape_id));

			// Read and set id
			buffer.unpack(p_id);
			particle->set_id(p_id);

			// Read all other quantities
			buffer.unpack(particle->get_position());
			buffer.unpack(particle->get_velocity());
			buffer.unpack(particle->get_angular_velocity());
			buffer.unpack(particle->get_orientation());
			buffer.unpack(particle->get_angular_velocity());
			buffer.unpack(particle->get_density());
			buffer.unpack(particle->get_force());
			buffer.unpack(particle->get_torque());
			buffer.unpack(particle->coll_force);
			buffer.unpack(particle->coll_force_last);
			buffer.unpack(particle->coll_torque);
			buffer.unpack(particle->coll_torque_last);

			particle->update();

			// Check that the file does intersect with the current domain
			if( ! geo::does_intersect(particle->get_bounding_sphere(), domain_bounds, arithmetic))
				continue;

			// Set processor id
			set_particle_proc_id(*particle);

			// Save the particle
			particles[p_id] = particle;
			add_solid_nodes_from_particle(*particle);
		}
	}

	is_init = true;
}

} /* namespace fsi */

} /* namespace plb */

#endif /* IMMERSEDBOUNDARYDYNAMICS_HH_ */
