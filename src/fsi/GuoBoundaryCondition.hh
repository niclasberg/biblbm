/*
 * GuoBoundaryCondition.hh
 *
 *  Created on: Aug 18, 2015
 *      Author: niber
 */

#ifndef GUOBOUNDARYCONDITION_HH_
#define GUOBOUNDARYCONDITION_HH_
#include "GuoBoundaryCondition.h"
#include "Boundary.h"

namespace plb {

namespace fsi {

template<class T, template<typename U> class Descriptor>
GuoRigidWallBoundary<T, Descriptor>::GuoRigidWallBoundary(const Boundary<T> & boundary, MultiBlockLattice3D<T, Descriptor> & lattice)
: boundary_(boundary), lattice_(lattice)
{
	applyProcessingFunctional(new GuoRigidWallBoundaryInstantiator<T, Descriptor>(*this), lattice_.getBoundingBox(), lattice_);
}

template<class T, template<typename U> class Descriptor>
void GuoRigidWallBoundary<T, Descriptor>::insert()
{
	integrateProcessingFunctional(new GuoRigidWallBoundaryFunctional<T, Descriptor>(*this), lattice_.getBoundingBox(), lattice_, 1);
}

template<class T, template<typename U> class Descriptor>
void GuoRigidWallBoundary<T, Descriptor>::put_node(Dot3D pos, Dot3D dir, T dist, const Array<T, 3> & wall_pos)
{
	Array<T, 3> dir_arr((T)dir.x, (T)dir.y, (T)dir.z);
	dir_arr /= norm(dir_arr);

	NodeInfo info;
	info.location = pos;
	info.delta = dist;
	info.direction = dir;
	info.wall_position = wall_pos;
	//std::cout << info.wall_position[0] << ", " << info.wall_position[1] << ", " << info.wall_position[2] << std::endl;
	info.velocity.resetToZero();
	wall_nodes_.push_back(info);
}

template<class T, template<typename U> class Descriptor>
void GuoRigidWallBoundary<T, Descriptor>::set_wall_velocity(const Array<T, 3> & velocity)
{
	for(int i = 0; i < wall_nodes_.size(); ++i)
		wall_nodes_[i].velocity = velocity;
}

template<class T, template<typename U> class Descriptor>
	template<class VelocityFunction> 
void GuoRigidWallBoundary<T, Descriptor>::set_wall_velocity(VelocityFunction velocity_function)
{
	for(int i = 0; i < wall_nodes_.size(); ++i)
		wall_nodes_[i].velocity = velocity_function(wall_nodes_[i].wall_position);
}

template<class T, template<typename U> class Descriptor>
void GuoRigidWallBoundaryInstantiator<T, Descriptor>::process(Box3D domain, BlockLattice3D<T, Descriptor> & lattice)
{
	Dot3D offset = lattice.getLocation();

	// Pre-compute the norm of the lattice vectors
	//Array<Array<T, 3>, Descriptor<T>::q> dirs;
	Array<T, Descriptor<T>::q> dir_norms;
	for(plint iDir = 1; iDir < Descriptor<T>::q; ++iDir) {
		dir_norms[iDir] = std::sqrt(Descriptor<T>::cNormSqr[iDir]);
	}

	// Locate wall nodes (i.e. nodes that are outside of the domain, but 
	// has a lattice vector that crosses the boundary)
	std::vector<Dot3D> wall_nodes;
	for(plint iX = domain.x0; iX <= domain.x1; ++iX) {
		for(plint iY = domain.y0; iY <= domain.y1; ++iY) {
			for(plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
				// Verify that the node is outside of the domain
				Array<T, 3> pos(iX+offset.x, iY+offset.y, iZ+offset.z);
				if( ! model_.boundary().contains(pos)) {
					// Check if any lattice vector crosses the surface of the boundary
					// if so, pick the direction of the crossing that best aligns with the normal of the surface.
					T cos_angle_min = 2;
					Dot3D dir_min(0, 0, 0);
					T delta_min = -1;
					Array<T, 3> wall_pos_min;
					wall_pos_min.resetToZero();
					bool success = false;

					for(plint iDir = 1; iDir < Descriptor<T>::q; ++iDir) {
						// Lattice vector
						Dot3D dir(Descriptor<T>::c[iDir][0], Descriptor<T>::c[iDir][1], Descriptor<T>::c[iDir][2]);
						Array<T, 3> dir2(dir.x, dir.y, dir.z);

						// Check if neighbour is a fluid node
						T t = -1;
						if(model_.boundary().trace_ray(pos, pos+dir2, t)) {
							//Normalized lattice vector
							dir2 /= dir_norms[iDir];

							//T dist = t * dir_norms[iDir]; //model_.boundary().distance_to_boundary(pos, dir2);
							Array<T, 3> wall_pos = pos + dir2*t;

							//if(dist >= 0 && dist <= dir_norms[iDir]) {
							Array<T, 3> normal = model_.boundary().get_normal(wall_pos);
							T cos_angle = dot(normal, dir2);
							if(cos_angle < cos_angle_min) {
								dir_min = dir;
								cos_angle_min = cos_angle;
								delta_min = 1 - t;//dist/dir_norms[iDir];
								wall_pos_min = wall_pos;
							}
							//}

							success = true;
						}
					}

					if(success) {
						//std::cout << cos_angle_max << ", " << delta_min << std::endl;
						//std::cout << pos[0] << ", " << pos[1] << ", " << pos[2] << " : "; 
						model_.put_node(Dot3D(iX+offset.x, iY+offset.y, iZ+offset.z), dir_min, delta_min, wall_pos_min);
					}
				}
			}
		}
	}

	/*std::fstream o("outside_pts.txt", std::ios::out);
	for(typename GuoRigidWallBoundary<T, Descriptor>::iterator it = model_.begin(); it != model_.end(); ++it)
		o << it->location.x << " " << it->location.y << " " << it->location.z
		  << " " << it->direction.x << " " << it->direction.y << " " << it->direction.z << std::endl;
	o.close();*/
}

template<class T, template<typename U> class Descriptor>
void GuoRigidWallBoundaryFunctional<T, Descriptor>::process(Box3D domain, BlockLattice3D<T, Descriptor> & lattice)
{
	Dot3D offset = lattice.getLocation();
	T rhoBar1, rhoBar2, rhoBar_w;
	Array<T, 3> j1, j2, j_w;
	Array<T, Descriptor<T>::q> fNeq1;
	Array<T, Descriptor<T>::q> fNeq2;
	Array<T, Descriptor<T>::q> fNeq_w;
	plint depth = 2;

	for(typename GuoRigidWallBoundary<T, Descriptor>::iterator it = model_.begin(); it != model_.end(); ++it) {
		Dot3D pos_rel(it->location.x - offset.x, it->location.y - offset.y, it->location.z - offset.z);

		plb::Cell<T, Descriptor> & cell = lattice.get(pos_rel.x, pos_rel.y, pos_rel.z);
		const plb::Cell<T, Descriptor> & cell1 = lattice.get(pos_rel.x + it->direction.x,
															 pos_rel.y + it->direction.y,
															 pos_rel.z + it->direction.z);
		const plb::Cell<T, Descriptor> & cell2 = lattice.get(pos_rel.x + 2*it->direction.x,
															 pos_rel.y + 2*it->direction.y,
															 pos_rel.z + 2*it->direction.z);

		// Correct the wall velocity (remove the force component)
		//for (int iD=0; iD<Descriptor<T>::d; ++iD)
		//	wall_vel[iD] = it->velocity[iD] - (T)0.5*getExternalForceComponent(cell,iD);

		// Evaluate non-equilibrium part of the distribution function
		cell1.getDynamics().computeRhoBarJ(cell1, rhoBar1, j1);
		cell2.getDynamics().computeRhoBarJ(cell2, rhoBar2, j2);
		T j1sqr = normSqr(j1);
		T j2sqr = normSqr(j2);
		for(plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
			fNeq1[iPop] = cell1[iPop] - cell1.getDynamics().computeEquilibrium(iPop, rhoBar1, j1, j1sqr);
			fNeq2[iPop] = cell2[iPop] - cell2.getDynamics().computeEquilibrium(iPop, rhoBar2, j2, j2sqr);
		}

		// Extrapolate wall momentum flux and non-equilibrium part of the distribution
		rhoBar_w = rhoBar1;
		Array<T,3> wall_j = Descriptor<T>::fullRho(rhoBar_w) * it->velocity;
		if (depth < 2) {
			j_w = (it->delta < (T)0.25) ? wall_j : (T)1./it->delta * (wall_j+(it->delta-(T)1.)*j1);
			fNeq_w = fNeq1;
		} else {  // depth >= 2
			if(it->delta < (T)0.75) {
				j_w = wall_j + (it->delta-(T)1.)*j1 + ((T)1.-it->delta)/((T)1.+it->delta)*((T)2.*wall_j+(it->delta-(T)1.)*j2);
				fNeq_w = it->delta*fNeq1 + ((T)1. - it->delta) * fNeq2;
			} else {
				j_w = (T)1./it->delta * (wall_j+(it->delta-(T)1.)*j1);
				fNeq_w = fNeq1;
			}
		}

		T jSqr_w = normSqr(j_w);
		for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
			cell[iPop] = cell.computeEquilibrium(iPop, rhoBar_w, j_w, jSqr_w)+fNeq_w[iPop];
		}

	}
}


template<class T, template<typename U> class Descriptor>
void GuoRigidWallBoundaryDebugger<T, Descriptor>::process(Box3D domain, ScalarField3D<T> & field)
{
	Dot3D offset = field.getLocation();
	for(plint iX = domain.x0; iX <= domain.x1; ++iX)
		for(plint iY = domain.y0; iY <= domain.y1; ++iY)
			for(plint iZ = domain.z0; iZ <= domain.z1; ++iZ)
				field.get(iX, iY, iZ) = 0;

	for(typename GuoRigidWallBoundary<T, Descriptor>::iterator it = model_.begin(); it != model_.end(); ++it) {
		Dot3D pos_rel(it->location.x - offset.x, it->location.y - offset.y, it->location.z - offset.z);

		field.get(pos_rel.x, pos_rel.y, pos_rel.z) = 1.;
		/*const plb::Cell<T, Descriptor> & cell1 = lattice.get(pos_rel.x + it->direction.x,
															 pos_rel.y + it->direction.y,
															 pos_rel.z + it->direction.z);
		const plb::Cell<T, Descriptor> & cell2 = lattice.get(pos_rel.x + 2*it->direction.x,
															 pos_rel.y + 2*it->direction.y,
															 pos_rel.z + 2*it->direction.z);*/
	}
}

}

}



#endif /* GUOBOUNDARYCONDITION_HH_ */
