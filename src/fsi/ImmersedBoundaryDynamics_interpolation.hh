#ifndef IMMERSEDBOUNDARYDYNAMICS_INTERPOLATION_HH_
#define IMMERSEDBOUNDARYDYNAMICS_INTERPOLATION_HH_
#include "ImmersedBoundaryDynamics.h"

namespace plb {

namespace fsi {

namespace interpolation {

namespace {

template<class T, class Arithmetic, class NodeType, class Dirac, class Interpolator>
void interpolate_impl(
		const Box3D & dom,
		TensorField3D<T, 3> & velocity,
		std::vector<NodeType> & node_container,
		const Arithmetic & arithmetic)
{
	Dot3D offset = velocity.getLocation();
	Box3D d_bounds;

	for(plint it = 0; it < node_container.size(); ++it) {
		// Dereference NodeType as a reference (it can be either a pointer or a reference)
		typename DeduceType<NodeType>::type & node = deref_maybe(node_container[it]);

		const Array<T, 3> pos_rel(node.pos[0] - offset.x, node.pos[1] - offset.y, node.pos[2] - offset.z);

		// Find compact support region of the Dirac function
		get_dirac_compact_support_box<T, Dirac>(pos_rel, d_bounds);

		// Interpolate fluid velocity
		node.vel.resetToZero();

		// The actual interpolation is outsourced to the Interpolator
		Interpolator::interpolate_node(d_bounds, dom, node, velocity, arithmetic, offset, pos_rel);
	}
}

// GeneralInterpolator:
// This interpolator will not try to interpolate outside of the domain specified by dom.
// It is thus useful for nodes that lie on the edge between two processor domains
template<class T, class Dirac>
struct GeneralInterpolator {
	template<class Arithmetic, class NodeType>
	static void interpolate_node(
			const Box3D & d_bounds,
			const Box3D & dom,
			NodeType & node,
			TensorField3D<T, 3> & velocity,
			const Arithmetic & arithmetic,
			const Dot3D & offset,
			const Array<T, 3> & pos_rel)
	{
		for(plint i = d_bounds.x0; i <= d_bounds.x1; ++i) {
			plint i2 = arithmetic.remap_index_x(i, offset.x);
			if(i2 < dom.x0 || i2 > dom.x1) continue;

			//T dirac_x = sampled_dirac.eval((T)i - pos_rel[0]);
			const T dirac_x = Dirac::eval((T)i - pos_rel[0]);
			for(plint j = d_bounds.y0; j <= d_bounds.y1; ++j) {
				plint j2 = arithmetic.remap_index_y(j, offset.y);
				if(j2 < dom.y0 || j2 > dom.y1) continue;

				const T dirac_xy = dirac_x * Dirac::eval((T)j - pos_rel[1]);
				//T dirac_xy = dirac_x * sampled_dirac.eval((T)j - pos_rel[1]);
				for(plint k = d_bounds.z0; k <= d_bounds.z1; ++k) {
					plint k2 = arithmetic.remap_index_z(k, offset.z);
					if(k2 < dom.z0 || k2 > dom.z1) continue;

					//T dirac_xyz = dirac_xy * sampled_dirac.eval((T)k - pos_rel[2]);
					const T dirac_xyz = dirac_xy * Dirac::eval((T)k - pos_rel[2]);
					node.vel += (dirac_xyz * velocity.get(i2, j2, k2));
				}
			}
		}
	}
};

// BulkInterpolator:
// This is an optimized version of the GeneralInterpolator, with the out-of-domain tests removed.
// This gives a much smoother loop for the processor without branching
template<class T, class Dirac>
struct BulkInterpolator {
	template<class Arithmetic, class NodeType>
	static void interpolate_node(
			const Box3D & d_bounds,
			const Box3D & dom,
			NodeType & node,
			TensorField3D<T, 3> & velocity,
			const Arithmetic & arithmetic,
			const Dot3D & offset,
			const Array<T, 3> & pos_rel)
	{
		for(plint i = d_bounds.x0; i <= d_bounds.x1; ++i) {
			const T dirac_x = Dirac::eval((T)i - pos_rel[0]);
			//T dirac_x = sampled_dirac.eval((T)i - pos_rel[0]);
			for(plint j = d_bounds.y0; j <= d_bounds.y1; ++j) {
				const T dirac_xy = dirac_x * Dirac::eval((T)j - pos_rel[1]);
				//T dirac_xy = dirac_x * sampled_dirac.eval((T)j - pos_rel[1]);
				for(plint k = d_bounds.z0; k <= d_bounds.z1; ++k) {
					const T dirac_xyz = dirac_xy * Dirac::eval((T)k - pos_rel[2]);
					//T dirac_xyz = dirac_xy * sampled_dirac.eval((T)k - pos_rel[2]);
					node.vel += (dirac_xyz * velocity.get(i, j, k));
				}
			}
		}
	}
};

// Optimized BulkInterpolator for RomaDirac
template<class T>
struct BulkInterpolator<T, RomaDirac<T, 3> > {
	template<class Arithmetic, class NodeType>
	static void interpolate_node(
			const Box3D & d_bounds,
			const Box3D & dom,
			NodeType & node,
			TensorField3D<T, 3> & velocity,
			const Arithmetic & arithmetic,
			const Dot3D & offset,
			const Array<T, 3> & pos_rel)
	{
		// Evaluate dirac values
		const T dx0 = RomaDirac<T, 3>::eval((T)d_bounds.x0 - pos_rel[0]);
		const T dx1 = RomaDirac<T, 3>::eval((T)d_bounds.x0 + (T)1. - pos_rel[0]);
		const T dx2 = RomaDirac<T, 3>::eval((T)d_bounds.x1 - pos_rel[0]);
		const T dy0 = RomaDirac<T, 3>::eval((T)d_bounds.y0 - pos_rel[1]);
		const T dy1 = RomaDirac<T, 3>::eval((T)d_bounds.y0 + (T)1. - pos_rel[1]);
		const T dy2 = RomaDirac<T, 3>::eval((T)d_bounds.y1 - pos_rel[1]);
		const T dz0 = RomaDirac<T, 3>::eval((T)d_bounds.z0 - pos_rel[2]);
		const T dz1 = RomaDirac<T, 3>::eval((T)d_bounds.z0 + (T)1. - pos_rel[2]);
		const T dz2 = RomaDirac<T, 3>::eval((T)d_bounds.z1 - pos_rel[2]);

		// Get indices
		// It will henceforth be assumed that d_bounds does not cross over a periodic edge,
		// then the indices to sample at are [i0, i0+1, i0+2] in x, and similarly for y and z.
		plint i0 = arithmetic.remap_index_x(d_bounds.x0, offset.x);
		plint j0 = arithmetic.remap_index_y(d_bounds.y0, offset.y);
		plint k0 = arithmetic.remap_index_z(d_bounds.z0, offset.z);;

		// Compute the interpolated value
		node.vel = dx0*(dy0*(dz0 * velocity.get(i0,   j0,   k0) + dz1 * velocity.get(i0,   j0,   k0+1) + dz2 * velocity.get(i0,   j0,   k0+2)) +
						dy1*(dz0 * velocity.get(i0,   j0+1, k0) + dz1 * velocity.get(i0,   j0+1, k0+1) + dz2 * velocity.get(i0,   j0+1, k0+2)) +
						dy2*(dz0 * velocity.get(i0,   j0+2, k0) + dz1 * velocity.get(i0,   j0+2, k0+1) + dz2 * velocity.get(i0,   j0+2, k0+2))) +
				   dx1*(dy0*(dz0 * velocity.get(i0+1, j0,   k0) + dz1 * velocity.get(i0+1, j0,   k0+1) + dz2 * velocity.get(i0+1, j0,   k0+2)) +
						dy1*(dz0 * velocity.get(i0+1, j0+1, k0) + dz1 * velocity.get(i0+1, j0+1, k0+1) + dz2 * velocity.get(i0+1, j0+1, k0+2)) +
						dy2*(dz0 * velocity.get(i0+1, j0+2, k0) + dz1 * velocity.get(i0+1, j0+2, k0+1) + dz2 * velocity.get(i0+1, j0+2, k0+2))) +
				   dx2*(dy0*(dz0 * velocity.get(i0+2, j0,   k0) + dz1 * velocity.get(i0+2, j0,   k0+1) + dz2 * velocity.get(i0+2, j0,   k0+2)) +
						dy1*(dz0 * velocity.get(i0+2, j0+1, k0) + dz1 * velocity.get(i0+2, j0+1, k0+1) + dz2 * velocity.get(i0+2, j0+1, k0+2)) +
						dy2*(dz0 * velocity.get(i0+2, j0+2, k0) + dz1 * velocity.get(i0+2, j0+2, k0+1) + dz2 * velocity.get(i0+2, j0+2, k0+2)));
	}
};

} /* namespace */

template<class T, class Arithmetic, class NodeType, class Dirac>
void interpolate_bulk(
		const Box3D & dom,
		TensorField3D<T, 3> & velocity,
		std::vector<NodeType> & node_container,
		const Arithmetic & arithmetic)
{
	interpolate_impl<T, Arithmetic, NodeType, Dirac, BulkInterpolator<T, Dirac> >(dom, velocity, node_container, arithmetic);
}

template<class T, class Arithmetic, class NodeType, class Dirac>
void interpolate(
		const Box3D & dom,
		TensorField3D<T, 3> & velocity,
		std::vector<NodeType> & node_container,
		const Arithmetic & arithmetic)
{
	interpolate_impl<T, Arithmetic, NodeType, Dirac, GeneralInterpolator<T, Dirac> >(dom, velocity, node_container, arithmetic);
}


// Interpolation for nodes close to the boundary
// We switch smoothly between the supplied Dirac function and a 2-point Dirac function
// when the distance to the boundary becomes less than the half support of the supplied Dirac
template<class T, class Arithmetic, class NodeType, class Dirac>
void interpolate_near_boundary(
		const Box3D & dom,
		TensorField3D<T, 3> & velocity,
		std::vector<NodeType> & node_container,
		const Arithmetic & arithmetic,
		Boundary<T> * boundary)
{
	Dot3D offset = velocity.getLocation();

	if(! boundary) {
		// No boundary, revert to the GeneralInterpolator
		interpolate<T, Arithmetic, NodeType, Dirac>(dom, velocity, node_container, arithmetic);
	} else {
		for(plint it = 0; it < node_container.size(); ++it) {
			// Dereference NodeType as a reference (it can be either a pointer or a reference)
			typename DeduceType<NodeType>::type & node = deref_maybe(node_container[it]);

			const Array<T, 3> pos_rel(node.pos[0] - offset.x, node.pos[1] - offset.y, node.pos[2] - offset.z);

			// Find compact support region of the Dirac function
			Box3D d_bounds, d_bounds_tophat;
			get_dirac_compact_support_box<T, Dirac>(pos_rel, d_bounds);
			get_dirac_compact_support_box<T, TopHatDirac<T, 3> >(pos_rel, d_bounds_tophat);

			// Compute distance to boundary and interpolation weight (linearly decreasing from 1 to 0)
			T weight = (boundary->distance_to_boundary(node.pos) - Dirac::half_support);

			if(weight >= 1) {
				// The node is far away from the boundary, use the supplied Dirac function
				node.vel.resetToZero();
				GeneralInterpolator<T, Dirac>::interpolate_node(d_bounds, dom, node, velocity, arithmetic, offset, pos_rel);
			} else if(weight <= 0) {
				// The node is close to the boundary, use the smallest possible Dirac function (top hat)
				node.vel.resetToZero();
				GeneralInterpolator<T, TopHatDirac<T, 3> >::interpolate_node(d_bounds_tophat, dom, node, velocity, arithmetic, offset, pos_rel);
			} else {
				// Compute contribution from the supplied Dirac
				node.vel.resetToZero();
				GeneralInterpolator<T, Dirac>::interpolate_node(d_bounds, dom, node, velocity, arithmetic, offset, pos_rel);
				Array<T, 3> suppliedDiracVelocity = node.vel;

				// Contribution from the TopHat dirac
				node.vel.resetToZero(); 
				GeneralInterpolator<T, TopHatDirac<T, 3> >::interpolate_node(d_bounds_tophat, dom, node, velocity, arithmetic, offset, pos_rel);
				node.vel = (1.-weight) * node.vel + weight*suppliedDiracVelocity;
			}
		}
	}
}

} /* namespace interpolation */

/********** Velocity interpolation **********/
template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::interpolate_velocity(
	const Box3D & dom,
	TensorField3D<T, 3> & velocity)
{
	Profile::start_timer("interpolate");

	// Interpolate non-local node velocity
	interpolation::interpolate_bulk<T, ArithmeticType, NonLocalNode<T>, Dirac>(
		dom, velocity, nonlocal_nodes, arithmetic);
	interpolation::interpolate<T, ArithmeticType, NonLocalNode<T>, Dirac>(
		dom, velocity, nonlocal_nodes_envelope, arithmetic);
	interpolation::interpolate_near_boundary<T, ArithmeticType, NonLocalNode<T>, Dirac>(
		dom, velocity, nonlocal_nodes_boundary, arithmetic, boundary);

	for(plint i = dom.x0; i <= dom.x1; ++i)
		for(plint j = dom.y0; j <= dom.y1; ++j)
			for(plint k = dom.z0; k <= dom.z1; ++k) {
				if(velocity.get(i, j, k)[0] < -1000)
					std::cout << "ASDFASDF";
			}

	// Send to appropriate processors
	send_interpolation_data();

	// Reset the velocity of all particle nodes to zero
	for(ObjMapIterator it = particles.begin(); it != particles.end(); ++it)
		for(plint i = 0; i < it->second->count_nodes(); ++i)
			it->second->get_node(i).vel.resetToZero();

	// Interpolate local
	interpolation::interpolate_bulk<T, ArithmeticType, Vertex<T> *, Dirac>(dom, velocity, local_nodes, arithmetic);
	interpolation::interpolate<T, ArithmeticType, Vertex<T> *, Dirac>(dom, velocity, local_nodes_envelope, arithmetic);
	interpolation::interpolate_near_boundary<T, ArithmeticType, Vertex<T> *, Dirac>(
			dom, velocity, local_nodes_boundary, arithmetic, boundary);
	Profile::stop_timer("interpolate");

	// Receive and reduce velocity
	Profile::start_timer("interpolate_reduction");
	receive_and_reduce_interpolation_data();
	Profile::stop_timer("interpolate_reduction");
}

namespace {
template<class T>
void pack_nodes(CommunicationBuffer & comm_buffer, std::vector<NonLocalNode<T> > & nodes)
{
	for(plint i = 0; i < nodes.size(); ++i) {
		const NonLocalNode<T> & node = nodes[i];
		//std::cout << "Sending " << node.particle_id << ", " << node.node_id << std::endl;
		comm_buffer.pack(node.proc_id, node.particle_id);
		comm_buffer.pack(node.proc_id, node.node_id);
		comm_buffer.pack(node.proc_id, node.vel);
	}
}
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::send_interpolation_data()
{
	comm_buffer.clear_send_buffer();
	pack_nodes(comm_buffer, nonlocal_nodes);
	pack_nodes(comm_buffer, nonlocal_nodes_envelope);
	pack_nodes(comm_buffer, nonlocal_nodes_boundary);
	comm_buffer.send_and_receive_no_wait(true);
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::receive_and_reduce_interpolation_data()
{
	// Synchronize and finalize data transfer
	comm_buffer.finalize_send_and_receive();

	// Perform a sum reduction of the velocity for each node
	plint p_id, node_id;
	Array<T, 3> vel;
	char * it = comm_buffer.recv_buffer_begin();

	while(it != comm_buffer.recv_buffer_end()) {
		comm_buffer.unpack(it, p_id);
		comm_buffer.unpack(it, node_id);
		comm_buffer.unpack(it, vel);

		ObjMapIterator p_it = particles.find(p_id);
		//std::cout << global::mpi().getRank() << ": PID " <<  p_id << " NID " << node_id << " vel (" << vel[0] << ", " << vel[1] << ", " << vel[2] << ")" << std::endl;
		if(p_it == particles.end())
			std::cerr << "ERROR: particle id = " << p_id << std::endl;

		PLB_PRECONDITION(p_it != particles.end())
		p_it->second->get_node(node_id).vel += vel;
	}
}
}

}



#endif /* IMMERSEDBOUNDARYDYNAMICS_INTERPOLATION_HH_ */
