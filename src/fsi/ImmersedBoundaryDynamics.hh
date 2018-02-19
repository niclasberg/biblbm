#ifndef IMMERSEDBOUNDARYDYNAMICS_HH_
#define IMMERSEDBOUNDARYDYNAMICS_HH_
#include "ImmersedBoundaryDynamics.h"
#include "atomicBlock/dataProcessor3D.h"
#include "core/geometry3D.h"
#include <algorithm>
#include "ParticleShapeFactory.h"
#include <sstream>
#include "MacroProcessors.h"
#include "Boundary.h"
#include <stdexcept>
#include "IO.h"
#include "TypeDeduction.h"

#define FSI_PROFILE
#include "Profile.h"


namespace plb {

namespace fsi {


/*------ ImmersedBoundaryDynamics3D ------*/
template<class T, template<typename U> class Descriptor, class Periodicity>
ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::ImmersedBoundaryDynamics3D(
		const MultiBlockLattice3D<T, Descriptor> & lattice,
		const ParticleShapeLibrary<T> & shape_library_)
: shape_library(shape_library_),
  management(lattice.getMultiBlockManagement()),
  num_particles(0),
  num_nodes(0),
  arithmetic(Periodicity::create_arithmetic(lattice.getBoundingBox())),
  global_bounding_box(lattice.getBoundingBox()),
  nx(lattice.getNx()),
  ny(lattice.getNy()),
  nz(lattice.getNz()),
  is_init_(false),
  z_buffer(arithmetic),
  boundary(0)
{
	initialize_domains();
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
	bulk_domain = domain;
	bulk_domain.enlarge_inplace(0.5 - Dirac::half_support - std::numeric_limits<T>::epsilon());

	std::vector<plint> proc_list;

	const SparseBlockStructure3D & block_structure = management.getSparseBlockStructure();
	const ThreadAttribution & thread_info = management.getThreadAttribution();

	// Find neighboring blocks
	std::vector<plint> neighboring_blocks;
	find_blocks_in_domain(domain_with_fsi_envelope, neighboring_blocks);

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
		box_with_envelope.enlarge_inplace(Dirac::half_support);

		MpiCommInfo<T> info;
		info.domain_with_fsi_envelope = box_with_envelope;
		info.domain = box;
		info.proc_id = proc_id;
		info.num_nodes_to_receive = 0;

		// Avoid adding duplicates
		if(std::find(proc_list.begin(), proc_list.end(), proc_id) == proc_list.end()) {
			proc_list.push_back(proc_id);
			comm_info[proc_id] = info;
		}
	}

	comm_buffer.set_proc_list(proc_list);
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

/* Destructor */
template<class T, template<typename U> class Descriptor, class Periodicity>
ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::~ImmersedBoundaryDynamics3D()
{
	clear();
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::clear()
{
	// Clear particles
	for(typename ObjMapType::iterator it = particles.begin(); it != particles.end(); ++it)
		delete it->second;
	particles.clear();
	num_particles = 0;
	particle_start_id_.clear();
	num_nodes = 0;

	// Clear nodes
	clear_nodes();
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::clear_nodes()
{
	local_nodes.clear();
	local_nodes_envelope.clear();
	local_nodes_boundary.clear();
	nonlocal_nodes.clear();
	nonlocal_nodes_envelope.clear();
	nonlocal_nodes_boundary.clear();
}

/********** Initialization **********/
template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::init()
{
	synchronize_particle_states();
	is_init_ = true;
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::synchronize_particle_states()
{
	std::vector<plint> particles_to_remove;
	comm_buffer.clear_send_buffer();

	pack_handoff_particles(particles_to_remove);
	pack_local_nodes();

	// Send data
	comm_buffer.send_and_receive_no_wait(true);

	// Clear old node lists and fill from all local particles
	clear_nodes();
	create_nodes();

	// Remove particles no longer owned by this domain
	for(std::vector<plint>::iterator it = particles_to_remove.begin(); it != particles_to_remove.end(); ++it) {
		delete particles.at(*it);
		particles.erase(*it);
	}

	// Receive data
	comm_buffer.finalize_send_and_receive();

	// Unpack received data
	char * it = comm_buffer.recv_buffer_begin();
	while(it != comm_buffer.recv_buffer_end()) {
		unpack_handoff_particles(it);
		unpack_nonlocal_nodes(it);
	}
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::pack_handoff_particles(std::vector<plint> & particles_to_remove)
{
	// We currently do not know how many particles and nodes we will send.
	// Therefore we pack Lazy values into the buffer which will be filled in later.
	std::map<plint, CommunicationBuffer::LazyValue<plint> > sizes;
	for(typename std::map<plint, MpiCommInfo<T> >::iterator it = comm_info.begin(); it != comm_info.end(); ++it) {
		sizes[it->first] = CommunicationBuffer::LazyValue<plint>(0);
		comm_buffer.pack(it->first, sizes.at(it->first));
	}

	// Particle hand-offs (particles that were local, but that now have their center of masses
	// in another processor)
	for(ObjMapIterator it = particles.begin(); it != particles.end(); ++it) {
		ParticleType & p = *(it->second);
		if( ! is_master_of(p)) {
			comm_buffer.delegate_pack(get_proc_id(p), p);
			sizes.at(get_proc_id(p)).get() += 1;
			particles_to_remove.push_back(it->first);
		}
	}

	// Finalize the lazy values
	for(std::map<plint, CommunicationBuffer::LazyValue<plint> >::iterator it = sizes.begin();
			it != sizes.end(); ++it) {
		//std::cout << it->first << " should receive " << it->second.get() << " particles" << std::endl;
		it->second.finalize();
	}
}

template<class T, template<typename U> class Descriptor, class Periodicity>
bool ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::is_master_of(
		const ParticleType & part) const
{
	return get_proc_id(part) == global::mpi().getRank();
}

template<class T, template<typename U> class Descriptor, class Periodicity>
plint ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::get_proc_id(
		const ParticleType & p) const
{
	// Get position in the periodic domain
	Array<T, 3> pos = p.center_of_mass();
	arithmetic.remap_position(pos);
	plint p_id = management.getSparseBlockStructure().
					locate(std::floor(pos[0]+0.5),
						   std::floor(pos[1]+0.5),
						   std::floor(pos[2]+0.5));
	if(p_id < 0)
		std::cerr << "Particle out of range: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
	return p_id;
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::pack_local_nodes()
{
	// We currently do not know how many particles nodes we will send, but we
	// want to put the number of nodes that will be received at the beginning of the
	// memory buffer. We therefore pack Lazy values into the buffer which will be filled in later.
	std::vector<CommunicationBuffer::LazyValue<plint> > sizes;
	for(typename std::map<plint, MpiCommInfo<T> >::iterator it = comm_info.begin(); it != comm_info.end(); ++it) {
		sizes.push_back(CommunicationBuffer::LazyValue<plint>(0));
		comm_buffer.pack(it->first, sizes[sizes.size()-1]);
	}

	for(ObjMapIterator it = particles.begin(); it != particles.end(); ++it) {
		ParticleType & p = *(it->second);
		plint proc_id = get_proc_id(p);

		// Pack nodes in neighboring processors
		plint size_i = 0;
		for(typename std::map<plint, MpiCommInfo<T> >::iterator it2 = comm_info.begin();
				it2 != comm_info.end(); ++it2, ++size_i) {
			// Don't send the nodes of particles that are handed off
			if(proc_id != it2->first) {
				// Do a quick test to eliminate particles whose bounding box
				// is not intersecting with the domain
				if(geo::does_intersect(p.bounding_box(), it2->second.domain_with_fsi_envelope, arithmetic)) {
					for(plint i = 0; i < p.count_nodes(); ++i) {
						if(it2->second.domain_with_fsi_envelope.contains(p.get_node(i).pos, arithmetic)) {
							comm_buffer.pack(it2->first, proc_id);
							comm_buffer.pack(it2->first, p.get_id());
							comm_buffer.pack(it2->first, i);
							comm_buffer.pack(it2->first, p.get_node(i).pos);
							comm_buffer.pack(it2->first, p.get_node(i).vel);
							comm_buffer.pack(it2->first, p.get_node(i).force);

							sizes[size_i].get() += 1;
						}
					}
				}
			}
		}
	}

	// Finalize the lazy values
	plint i = 0;
	for(typename std::map<plint, MpiCommInfo<T> >::iterator it = comm_info.begin();
			it != comm_info.end(); ++it, ++i) {
		it->second.num_nodes_to_receive = sizes[i].get();
		sizes[i].finalize();
	}
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::create_nodes()
{
	const Array<T, 3> mid(0.5*(domain.x0 + domain.x1),
						  0.5*(domain.y0 + domain.y1),
						  0.5*(domain.z0 + domain.z1));
	for(ObjMapIterator it = particles.begin(); it != particles.end(); ++it) {
		ParticleType & p = *(it->second);
		plint proc_id = get_proc_id(p);

		if(proc_id == global::mpi().getRank()) {
			create_nodes_from_particle(p);
		} else {
			for(plint i = 0; i < p.count_nodes(); ++i) {
				Vertex<T> & n = p.get_node(i);
				Array<T, 3> pos = n.pos;
				arithmetic.shift_periodically_to_minimize_distance_to(mid, pos);

				if(domain_with_fsi_envelope.contains(pos, arithmetic)) {
					NonLocalNode<T> node;
					node.proc_id = proc_id;
					node.node_id = i;
					node.particle_id = p.get_id();
					node.pos = pos;
					node.vel = n.vel;
					node.force = n.force;
					add_nonlocal_node(node);
				}
			}
		}
	}
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::add_nonlocal_node(NonLocalNode<T> & node)
{
	if(boundary && boundary->distance_to_boundary_less_than(node.pos, Dirac::half_support))
		nonlocal_nodes_boundary.push_back(node);
	else if(bulk_domain.contains(node.pos, arithmetic))
		nonlocal_nodes.push_back(node);
	else
		nonlocal_nodes_envelope.push_back(node);
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::create_nodes_from_particle(ParticleType & p)
{
	for(typename ParticleBase3D<T>::vertex_iterator it = p.begin(); it != p.end(); ++it) {
		if(domain_with_fsi_envelope.contains(it->pos, arithmetic)) {
			if(boundary && boundary->distance_to_boundary_less_than(it->pos, Dirac::half_support))
				local_nodes_boundary.push_back(&(*it));
			else if(bulk_domain.contains(it->pos, arithmetic))
				local_nodes.push_back(&(*it));
			else
				local_nodes_envelope.push_back(&(*it));
		}
	}
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::unpack_handoff_particles(char *& it)
{
	// Read number of particles
	plint num_part;
	comm_buffer.unpack(it, num_part);

	// Unpack particles
	for(plint i = 0; i < num_part; ++i) {
		ParticleBase3D<T> * p = particleFactory<T>().create(it, shape_library);

		// Create node pointers
		create_nodes_from_particle(*p);
		particles[p->get_id()] = p;
	}
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::unpack_nonlocal_nodes(char *& it)
{
	const Array<T, 3> mid(0.5*(domain.x0 + domain.x1),
						  0.5*(domain.y0 + domain.y1),
						  0.5*(domain.z0 + domain.z1));

	// Read number of nodes
	plint num_nodes;
	comm_buffer.unpack(it, num_nodes);

	// Unpack nodes
	for(plint i = 0; i < num_nodes; ++i) {
		NonLocalNode<T> node;
		comm_buffer.unpack(it, node.proc_id);
		comm_buffer.unpack(it, node.particle_id);
		comm_buffer.unpack(it, node.node_id);
		comm_buffer.unpack(it, node.pos);
		comm_buffer.unpack(it, node.vel);
		comm_buffer.unpack(it, node.force);

		arithmetic.shift_periodically_to_minimize_distance_to(mid, node.pos);
		if(node.proc_id != global::mpi().getRank())
			add_nonlocal_node(node);
	}
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::init(
		const Box3D & dom,
		BlockLattice3D<T, Descriptor> & lattice)
{
	voxelize();
	synchronize_particle_states_and_voxelization(dom, lattice);
	apply_voxelization(dom, lattice);
	is_init_ = true;
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::synchronize_particle_states_and_voxelization(
		const Box3D & dom,
		BlockLattice3D<T, Descriptor> & lattice)
{
	std::vector<plint> particles_to_remove;
	comm_buffer.clear_send_buffer();

	pack_handoff_particles(particles_to_remove);
	pack_local_nodes();
	pack_voxelization_data();

	// Send and receive
	comm_buffer.send_and_receive_no_wait(true);

	// Clear old node lists and fill from all local particles
	clear_nodes();
	create_nodes();

	// Remove particles no longer owned by this domain
	for(std::vector<plint>::iterator it = particles_to_remove.begin(); it != particles_to_remove.end(); ++it) {
		delete particles.at(*it);
		particles.erase(*it);
	}

	// Receive data
	comm_buffer.finalize_send_and_receive();

	// Unpack received data
	char * it = comm_buffer.recv_buffer_begin();
	while(it != comm_buffer.recv_buffer_end()) {
		unpack_handoff_particles(it);
		unpack_nonlocal_nodes(it);
		unpack_voxelization_data(it);
	}
}

/******** Vertex motion ********/
template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::move_vertices()
{
	move_vertices_impl();
	synchronize_particle_states();
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::move_vertices_and_revoxelize(
		const Box3D & dom,
		BlockLattice3D<T, Descriptor> & lattice)
{
	Profile::start_timer("move");
	// Move local vertices
	move_vertices_impl();
	Profile::stop_timer("move");

	// Voxelize the particles
	Profile::start_timer("voxelizer_rehash");
	voxelize();
	Profile::stop_timer("voxelizer_rehash");

	// Synchronize particle states
	Profile::start_timer("synch");
	synchronize_particle_states_and_voxelization(dom, lattice);
	Profile::stop_timer("synch");

	Profile::start_timer("voxelizer_apply");
	apply_voxelization(dom, lattice);
	Profile::stop_timer("voxelizer_apply");
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::move_vertices_impl()
{
	// Integrate the equations of motion for all particles
	for(ObjMapIterator it = particles.begin(); it != particles.end(); ++it) {
		it->second->move_vertices(boundary);
		it->second->map_center_of_mass_to_periodic_grid(arithmetic);
	}
}

/*
 * add_particle
 * Adds a particle to the simulation. Particles not in the domain are discarded.
 */
template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::add_particle(const ParticleType * particle)
{
	// Create an id for the particle.
	plint p_id = get_next_id();
	add_particle(particle, p_id);
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::add_particle(const ParticleType * particle, plint p_id)
{
	// For now, just save the number of nodes in the start id array
	particle_start_id_[p_id] = num_nodes;
	num_nodes += particle->count_nodes();

	// Only add particles that intersects this domain
	if(is_master_of(*particle)) {
		// Copy the particle
		ParticleType * p_copy = particle->clone();

		// Set particle id
		p_copy->map_center_of_mass_to_periodic_grid(arithmetic);
		p_copy->set_id(p_id);

		particles[p_id] = p_copy;
	}
}

template<class T, template<typename U> class Descriptor, class Periodicity>
bool ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::get_particle(plint i, ParticleType *& particle)
{
	ObjMapIterator it = particles.find(i);
	if(it == particles.end())
		return false;
	particle = it->second;
	return true;
}

template<class T, template<typename U> class Descriptor, class Periodicity>
bool ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::get_particle(plint i, const ParticleType *& particle) const
{
	ObjMapConstIterator it = particles.find(i);
	if(it == particles.end())
		return false;
	particle = it->second;
	return true;
}

/*** Global reductions ***/
template<class T, template<typename U> class Descriptor, class Periodicity>
plint ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::count_total_particles() const
{
	plint local_particles = particles.size();
	plint total_particles;
	global::mpi().reduce(local_particles, total_particles, MPI_SUM);

	return total_particles;
}

template<class T, template<typename U> class Descriptor, class Periodicity>
T ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::compute_total_volume() const
{
	T local_volume = T();
	T total_volume = 0;;
	for(typename ObjMapType::const_iterator it = particles.begin(); it != particles.end(); ++it)
		local_volume += it->second->volume();

	global::mpi().reduce(local_volume, total_volume, MPI_SUM);
	return total_volume;
}


} /* namespace fsi */

} /* namespace plb */

#include "ImmersedBoundaryDynamics_interpolation.hh"
#include "ImmersedBoundaryDynamics_io.hh"
#include "ImmersedBoundaryDynamics_forces.hh"
#include "ImmersedBoundaryDynamics_voxelization.hh"

#endif /* IMMERSEDBOUNDARYDYNAMICS_HH_ */
