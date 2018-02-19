/*
 * ImmersedBoundaryDynamics_forces.hh
 *
 *  Created on: Apr 20, 2016
 *      Author: niber
 */

#ifndef IMMERSEDBOUNDARYDYNAMICS_FORCES_HH_
#define IMMERSEDBOUNDARYDYNAMICS_FORCES_HH_
#include "ImmersedBoundaryDynamics.h"
#include "CollisionForces.h"
#include <tr1/unordered_map>

namespace plb {

namespace fsi {

namespace spreading {

namespace {

template<class T, template<typename U> class Descriptor, class Dirac, class Arithmetic, class NodeType, class Spreader>
void spread_forces_impl(
		const Box3D & dom,
		BlockLattice3D<T, Descriptor> & lattice,
		std::vector<NodeType> & node_container,
		const Arithmetic & arithmetic,
		const SampledDirac<T, Dirac, 3> & sampled_dirac)
{
	Dot3D offset = lattice.getLocation();
	Box3D d_bounds;
	T dirac_x, dirac_xy, dirac_xyz;

	for(plint it = 0; it < node_container.size(); ++it) {
		typename DeduceType<NodeType>::type & node = deref_maybe(node_container[it]);

		const Array<T, 3> pos_rel(node.pos[0] - offset.x, node.pos[1] - offset.y, node.pos[2] - offset.z);

		// Find compact support region of the Dirac function
		get_dirac_compact_support_box<T, Dirac>(pos_rel, d_bounds);

		Spreader::spread(d_bounds, pos_rel, dom, lattice, node, arithmetic, sampled_dirac, offset);
	}
}

template<class T>
void add_to_cArray(T * p, const Array<T, 3> & array)
{
	p[0] += array[0];
	p[1] += array[1];
	p[2] += array[2];
}

template<class T, class Dirac>
struct GeneralSpreading {
	template<template<typename U> class Descriptor, class Arithmetic, class NodeType>
	static void spread(
			const Box3D & d_bounds,
			const Array<T, 3> & pos_rel,
			const Box3D & dom,
			BlockLattice3D<T, Descriptor> & lattice,
			NodeType & node,
			const Arithmetic & arithmetic,
			const SampledDirac<T, Dirac, 3> & sampled_dirac,
			const Dot3D & offset)
	{
		for(plint i = d_bounds.x0; i <= d_bounds.x1; ++i) {
			plint i2 = arithmetic.remap_index_x(i, offset.x);
			if(i2 < dom.x0 || i2 > dom.x1) continue;

			T dirac_x = sampled_dirac.eval((T)i - pos_rel[0]);
			for(plint j = d_bounds.y0; j <= d_bounds.y1; ++j) {
				plint j2 = arithmetic.remap_index_y(j, offset.y);
				if(j2 < dom.y0 || j2 > dom.y1) continue;

				T dirac_xy = dirac_x * sampled_dirac.eval((T)j - pos_rel[1]);
				for(plint k = d_bounds.z0; k <= d_bounds.z1; ++k) {
					plint k2 = arithmetic.remap_index_z(k, offset.z);
					if(k2 < dom.z0 || k2 > dom.z1) continue;
					T dirac_xyz = dirac_xy * sampled_dirac.eval((T)k - pos_rel[2]);

					add_to_cArray(lattice.get(i2, j2, k2).getExternal(Descriptor<T>::ExternalField::forceBeginsAt),
								dirac_xyz * node.force);
				}
			}
		}
	}
};

// Optimized spreading (no out of range checks)
template<class T, class Dirac>
struct BulkSpreading {
	template<template<typename U> class Descriptor, class Arithmetic, class NodeType>
	static void spread(
			const Box3D & d_bounds,
			const Array<T, 3> & pos_rel,
			const Box3D & dom,
			BlockLattice3D<T, Descriptor> & lattice,
			NodeType & node,
			const Arithmetic & arithmetic,
			const SampledDirac<T, Dirac, 3> & sampled_dirac,
			const Dot3D & offset)
	{
		for(plint i = d_bounds.x0; i <= d_bounds.x1; ++i) {
			plint i2 = arithmetic.remap_index_x(i, offset.x);
			T dirac_x = sampled_dirac.eval((T)i - pos_rel[0]);
			for(plint j = d_bounds.y0; j <= d_bounds.y1; ++j) {
				plint j2 = arithmetic.remap_index_y(j, offset.y);
				T dirac_xy = dirac_x * sampled_dirac.eval((T)j - pos_rel[1]);
				for(plint k = d_bounds.z0; k <= d_bounds.z1; ++k) {
					plint k2 = arithmetic.remap_index_z(k, offset.z);
					T dirac_xyz = dirac_xy * sampled_dirac.eval((T)k - pos_rel[2]);
					add_to_cArray(lattice.get(i2, j2, k2).getExternal(Descriptor<T>::ExternalField::forceBeginsAt),
								dirac_xyz * node.force);
				}
			}
		}
	}
};

// Optimized spreading for Roma Dirac
template<class T>
struct BulkSpreading<T, RomaDirac<T, 3> > {
	template<template<typename U> class Descriptor, class Arithmetic, class NodeType>
	static void spread(
			const Box3D & d_bounds,
			const Array<T, 3> & pos_rel,
			const Box3D & dom,
			BlockLattice3D<T, Descriptor> & lattice,
			NodeType & node,
			const Arithmetic & arithmetic,
			const SampledDirac<T, RomaDirac<T, 3>, 3> & sampled_dirac,
			const Dot3D & offset)
	{
		// Evaluate dirac values
		const T dx[3] = {sampled_dirac.eval((T)d_bounds.x0 - pos_rel[0]),
						 sampled_dirac.eval((T)d_bounds.x0 + (T)1. - pos_rel[0]),
						 sampled_dirac.eval((T)d_bounds.x1 - pos_rel[0])};
		const T dy0 = sampled_dirac.eval((T)d_bounds.y0 - pos_rel[1]);
		const T dy1 = sampled_dirac.eval((T)d_bounds.y0 + (T)1. - pos_rel[1]);
		const T dy2 = sampled_dirac.eval((T)d_bounds.y1 - pos_rel[1]);
		const T dz0 = sampled_dirac.eval((T)d_bounds.z0 - pos_rel[2]);
		const T dz1 = sampled_dirac.eval((T)d_bounds.z0 + (T)1. - pos_rel[2]);
		const T dz2 = sampled_dirac.eval((T)d_bounds.z1 - pos_rel[2]);

		// Get indices
		// It will henceforth be assumed that d_bounds does not cross over a periodic edge,
		// then the indices to spread to are [i0, i0+1, i0+2] in x, and similarly for y and z.
		plint i0 = arithmetic.remap_index_x(d_bounds.x0, offset.x);
		plint j0 = arithmetic.remap_index_y(d_bounds.y0, offset.y);
		plint k0 = arithmetic.remap_index_z(d_bounds.z0, offset.z);
		plint forceInd = Descriptor<T>::ExternalField::forceBeginsAt;

		// Spead force
		for(plint x = 0; x < 3; ++x) {
			// Unroll the two inner loops
			add_to_cArray(lattice.get(i0+x, j0,   k0  ).getExternal(forceInd), dx[x] * dy0 * dz0 * node.force);
			add_to_cArray(lattice.get(i0+x, j0,   k0+1).getExternal(forceInd), dx[x] * dy0 * dz1 * node.force);
			add_to_cArray(lattice.get(i0+x, j0,   k0+2).getExternal(forceInd), dx[x] * dy0 * dz2 * node.force);
			add_to_cArray(lattice.get(i0+x, j0+1, k0  ).getExternal(forceInd), dx[x] * dy1 * dz0 * node.force);
			add_to_cArray(lattice.get(i0+x, j0+1, k0+1).getExternal(forceInd), dx[x] * dy1 * dz1 * node.force);
			add_to_cArray(lattice.get(i0+x, j0+1, k0+2).getExternal(forceInd), dx[x] * dy1 * dz2 * node.force);
			add_to_cArray(lattice.get(i0+x, j0+2, k0  ).getExternal(forceInd), dx[x] * dy2 * dz0 * node.force);
			add_to_cArray(lattice.get(i0+x, j0+2, k0+1).getExternal(forceInd), dx[x] * dy2 * dz1 * node.force);
			add_to_cArray(lattice.get(i0+x, j0+2, k0+2).getExternal(forceInd), dx[x] * dy2 * dz2 * node.force);
		}
	}
};

} /* namespace */

template<class T, template<typename U> class Descriptor, class Dirac, class Arithmetic, class NodeType>
void spread_forces(
		const Box3D & dom,
		BlockLattice3D<T, Descriptor> & lattice,
		std::vector<NodeType> & node_container,
		const Arithmetic & arithmetic,
		const SampledDirac<T, Dirac, 3> & sampled_dirac)
{
	spread_forces_impl<T, Descriptor, Dirac, Arithmetic, NodeType, GeneralSpreading<T, Dirac> >(dom, lattice, node_container, arithmetic, sampled_dirac);
}

template<class T, template<typename U> class Descriptor, class Dirac, class Arithmetic, class NodeType>
void spread_forces_bulk(
		const Box3D & dom,
		BlockLattice3D<T, Descriptor> & lattice,
		std::vector<NodeType> & node_container,
		const Arithmetic & arithmetic,
		const SampledDirac<T, Dirac, 3> & sampled_dirac)
{
	spread_forces_impl<T, Descriptor, Dirac, Arithmetic, NodeType, BulkSpreading<T, Dirac> >(dom, lattice, node_container, arithmetic, sampled_dirac);
}

// Boundary spreading. Here we account for that some grid points may lie outside of the domain.
// A Dirac function will be formed including only the points that are within the domain.
// The remaining points will be weighted such that the moment conditions:
//   sum dirac(x-X) = 1
//   sum (x-X)*dirac(x-X) = 0
// still hold.
template<class T, template<typename U> class Descriptor, class Dirac, class Arithmetic, class NodeType>
void spread_forces_near_boundary(
		const Box3D & dom,
		BlockLattice3D<T, Descriptor> & lattice,
		std::vector<NodeType> & node_container,
		const Arithmetic & arithmetic,
		Boundary<T> * boundary)
{
	Dot3D offset = lattice.getLocation();
	Box3D d_bounds;

	for(plint it = 0; it < node_container.size(); ++it) {
		// Dereference NodeType as a reference (it can be either a pointer or a reference)
		typename DeduceType<NodeType>::type & node = deref_maybe(node_container[it]);

		// Find compact support region of the Dirac function
		get_dirac_compact_support_box<T, Dirac>(node.pos, d_bounds);

		// Determine node topology (find points that are outside of the domain)
		DiracWithMissingPoints<T, Dirac> dirac(node.pos);
		for(plint i = d_bounds.x0; i <= d_bounds.x1; ++i) {
			plint i2 = arithmetic.remap_index_x(i);
			for(plint j = d_bounds.y0; j <= d_bounds.y1; ++j) {
				plint j2 = arithmetic.remap_index_y(j);
				for(plint k = d_bounds.z0; k <= d_bounds.z1; ++k) {
					plint k2 = arithmetic.remap_index_z(k);
					if( ! boundary->contains(Array<T, 3>((T)i2, (T)j2, (T)k2))) {
						dirac.setNodeIsValid(i, j, k, false);
					}
				}
			}
		}

		// Construct Dirac function
		dirac.computeWeights();

		// Spread forces
		for(plint i=0; i < dirac.count_points(); ++i) {
			Dot3D p = dirac.get_dirac_point(i).node_pos;
			p -= offset;
			plint i2 = arithmetic.remap_index_x(p.x, offset.x);
			plint j2 = arithmetic.remap_index_y(p.y, offset.y);
			plint k2 = arithmetic.remap_index_z(p.z, offset.z);

			if(contained(i2, j2, k2, dom)) {
				add_to_cArray(lattice.get(i2, j2, k2).getExternal(Descriptor<T>::ExternalField::forceBeginsAt),
						dirac.get_dirac_point(i).weight * dirac.get_dirac_point(i).dirac_val * node.force);
			}
		}
	}
}

} /* namespace spreading */

/******** Compute and spread forces ********/
template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::compute_and_spread_forces(
		const Box3D & dom,
		BlockLattice3D<T, Descriptor> & lattice)
{
	// Compute forces
	Profile::start_timer("compute_forces");
	compute_forces();
	Profile::stop_timer("compute_forces");

	// Send non-local nodes
	Profile::start_timer("spread_forces");
	comm_buffer.clear_send_buffer();
	pack_local_particle_forces();
	if( ! pp_forces_.empty())
		pack_nonlocal_node_forces();
	comm_buffer.send_and_receive_no_wait(true);

	// Spread local force
	spreading::spread_forces_bulk(dom, lattice, local_nodes, arithmetic, sampled_dirac);

	// Receive non-local nodes
	comm_buffer.finalize_send_and_receive();
	for(char * it = comm_buffer.recv_buffer_begin(); it != comm_buffer.recv_buffer_end(); ) {
		unpack_local_particle_forces(it);
		if( ! pp_forces_.empty())
			reduce_nonlocal_node_forces(it);
	}

	// Spread non-local forces and forces in the envelope and boundary
	spreading::spread_forces(dom, lattice, local_nodes_envelope, arithmetic, sampled_dirac);
	spreading::spread_forces_bulk(dom, lattice, nonlocal_nodes, arithmetic, sampled_dirac);
	spreading::spread_forces(dom, lattice, nonlocal_nodes_envelope, arithmetic, sampled_dirac);
	spreading::spread_forces_near_boundary<T, Descriptor, Dirac, ArithmeticType, Vertex<T> * >(dom, lattice, local_nodes_boundary, arithmetic, boundary);
	spreading::spread_forces_near_boundary<T, Descriptor, Dirac, ArithmeticType, NonLocalNode<T> >(dom, lattice, nonlocal_nodes_boundary, arithmetic, boundary);
	Profile::stop_timer("spread_forces");
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::compute_forces()
{
	// Set all forces to zero
	set_forces_to_zero();

	// Compute forces for all particles (updates local nodes)
	for(ObjMapIterator it = particles.begin(); it != particles.end(); ++it) {
		for(int i = 0; i < force_decorators_.size(); ++i)
			force_decorators_[i]->apply_force(it->second);
		it->second->compute_forces();
	}

	// Compute particle-particle interaction forces for both local and non-local nodes
	for(int i = 0; i < pp_forces_.size(); ++i)
		compute_particle_particle_interaction_impl(*(pp_forces_[i]), 0);
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::set_forces_to_zero()
{
	// Reset local particle forces
	for(ObjMapIterator it = particles.begin(); it != particles.end(); ++it)
		it->second->reset_forces();

	// Reset forces on non-local nodes
	for(plint i = 0; i < nonlocal_nodes.size(); ++i)
		nonlocal_nodes[i].force.resetToZero();
	for(plint i = 0; i < nonlocal_nodes_envelope.size(); ++i)
		nonlocal_nodes_envelope[i].force.resetToZero();
	for(plint i = 0; i < nonlocal_nodes_boundary.size(); ++i)
		nonlocal_nodes_boundary[i].force.resetToZero();
}

namespace detail {
template<class T, class CollHandler>
void add_nonlocal_nodes_to_collision_handler(std::vector<NonLocalNode<T> > & nodes, CollHandler & coll_handler)
{
	for(plint i = 0; i < nodes.size(); ++i) {
		coll_handler.add_node(nodes[i].pos, &(nodes[i].force), nodes[i].particle_id);
	}
}
} /* namespace detail */

template<class T, template<typename U> class Descriptor, class Periodicity>
	template<class InteractionFunctional>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::compute_particle_particle_interaction_impl(
		const InteractionFunctional & interaction,
		const Boundary<T> * boundary)
{
	// Domain for collision checking
	geo::Rect<T> bb(
			domain.x0-0.5, domain.x1+0.5-std::numeric_limits<T>::epsilon(),
			domain.y0-0.5, domain.y1+0.5-std::numeric_limits<T>::epsilon(),
			domain.z0-0.5, domain.z1+0.5-std::numeric_limits<T>::epsilon());

	geo::Rect<T> coll_domain = bb.enlarge(interaction.get_cutoff_distance());

	// We only communicate nodes that are at a distance of the half_support of the Dirac function
	// to neighboring processors, the interaction calculation will be errounous if the cutoff distance
	// is larger than this.
	if(interaction.get_cutoff_distance() > (Dirac::half_support - 0.5)) {
		std::cerr << "The interaction distance cannot be greater than the lattice spacing (1)" << std::endl;
		exit(-1);
	}

	collision::CollisionHandler<T, typename Periodicity::ArithmeticType, InteractionFunctional>
		collision_handler(bb, arithmetic, interaction);

	// Add non-local nodes
	detail::add_nonlocal_nodes_to_collision_handler(nonlocal_nodes, collision_handler);
	detail::add_nonlocal_nodes_to_collision_handler(nonlocal_nodes_envelope, collision_handler);

	// Add local nodes
	for(ObjMapIterator it = particles.begin(); it != particles.end(); ++it) {
		ParticleBase3D<T> * p = it->second;
		for(plint i = 0; i < p->count_nodes(); ++i) {
			if(coll_domain.contains(p->get_node(i).pos))
				collision_handler.add_node(p->get_node(i).pos, &(p->get_node(i).force), p->get_id());
		}
	}

	// Evaluate collision forces on all nodes
	collision_handler.compute_collision_forces();

	// Compute wall collision forces if a boundary is supplied
	if(boundary)
		collision_handler.compute_wall_collision_forces(*boundary);
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::pack_local_particle_forces()
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
							comm_buffer.pack(it2->first, p.get_id());
							comm_buffer.pack(it2->first, i);
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
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::pack_nonlocal_node_forces()
{
	// Just save pointers to the nonlocal node arrays to avoid code duplication
	std::vector<NonLocalNode<T> > * containers[3] = {&nonlocal_nodes, &nonlocal_nodes_envelope, &nonlocal_nodes_boundary };

	// Count the number of nodes to send to each proc
	std::map<plint, plint> sizes;
	for(typename std::map<plint, MpiCommInfo<T> >::iterator it = comm_info.begin(); it != comm_info.end(); ++it)
		sizes[it->first] = 0;

	for(plint ci = 0; ci < 3; ++ci) {
		std::vector<NonLocalNode<T> > & nodes = *(containers[ci]);
		for(plint i = 0; i < nodes.size(); ++i) {
			sizes.at(nodes[i].proc_id) += 1;
		}
	}

	// Pack number of nodes to send
	for(std::map<plint, plint>::iterator it = sizes.begin(); it != sizes.end(); ++it)
		comm_buffer.pack(it->first, it->second);

	// Pack data
	for(plint ci = 0; ci < 3; ++ci) {
		std::vector<NonLocalNode<T> > & nodes = *(containers[ci]);
		for(plint i = 0; i < nodes.size(); ++i) {
			const NonLocalNode<T> & node = nodes[i];
			comm_buffer.pack(node.proc_id, node.particle_id);
			comm_buffer.pack(node.proc_id, node.node_id);
			comm_buffer.pack(node.proc_id, node.force);
		}
	}
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::unpack_local_particle_forces(char *& it)
{
	// Put pointers to the nodes in a hash map for faster lookup
	std::tr1::unordered_map<plint, NonLocalNode<T> *> nonlocal_node_map;
	std::vector<NonLocalNode<T> > * containers[3] = {&nonlocal_nodes, &nonlocal_nodes_envelope , &nonlocal_nodes_boundary };
	for(int ci = 0; ci < 3; ++ci) {
		std::vector<NonLocalNode<T> > & nodes = *(containers[ci]);
		for(plint i = 0; i < nodes.size(); ++i) {
			NonLocalNode<T> & node = nodes[i];
			nonlocal_node_map[particle_start_id_[node.particle_id] + node.node_id] = &node;
		}
	}

	plint num_nodes, pid, nid;
	Array<T, 3> force;
	comm_buffer.unpack(it, num_nodes);
	for(plint i = 0; i < num_nodes; ++i) {
		comm_buffer.unpack(it, pid);
		comm_buffer.unpack(it, nid);
		comm_buffer.unpack(it, force);
		nonlocal_node_map.find(particle_start_id_[pid] + nid)->second->force += force;
	}
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::reduce_nonlocal_node_forces(char *& it)
{
	plint p_id, node_id;
	Array<T, 3> force;
	plint num_nodes;
	comm_buffer.unpack(it, num_nodes);

	for(plint i = 0; i < num_nodes; ++i) {
		comm_buffer.unpack(it, p_id);
		comm_buffer.unpack(it, node_id);
		comm_buffer.unpack(it, force);

		ObjMapIterator p_it = particles.find(p_id);
		//std::cout << global::mpi().getRank() << ": PID " <<  p_id << " NID " << node_id << " vel (" << vel[0] << ", " << vel[1] << ", " << vel[2] << ")" << std::endl;
		if(p_it == particles.end())
			std::cerr << "ERROR: particle id = " << p_id << std::endl;
		p_it->second->get_node(node_id).force += force;
	}
}

/****** Collision forces ********/
template<class T, template<typename U> class Descriptor, class Periodicity>
	template<class InteractionFunctional>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::compute_collision_forces(
		const InteractionFunctional & interaction,
		const Boundary<T> * boundary)
{
	set_forces_to_zero();
	compute_particle_particle_interaction_impl(interaction, boundary);

	// Communicate forces on nonlocal nodes
	comm_buffer.clear_send_buffer();
	pack_nonlocal_node_forces();
	comm_buffer.send_and_receive_no_wait(true);
	comm_buffer.finalize_send_and_receive();

	// Reduce the total force on all local particles
	for(char * it = comm_buffer.recv_buffer_begin(); it != comm_buffer.recv_buffer_end(); )
		reduce_nonlocal_node_forces(it);
}

}
}


#endif /* IMMERSEDBOUNDARYDYNAMICS_FORCES_HH_ */
