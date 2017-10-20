/*
 * DataProcessors.hh
 *
 *  Created on: Jun 24, 2015
 *      Author: niber
 */

#ifndef DATAPROCESSORS_HH_
#define DATAPROCESSORS_HH_
#include "DataProcessors.h"

namespace plb {

namespace fsi {

/*** ComputeConcentrationFunctional ***/
template<class T, template<typename U> class Descriptor>
void ComputeConcentrationFunctional<T, Descriptor>::process(
		Box3D domain,
		BlockLattice3D<T, Descriptor> & lattice,
		ScalarField3D<T> & conc)
{
	Dot3D offset = computeRelativeDisplacement(lattice, conc);

	for(plint i = domain.x0; i <= domain.x1; ++i)
		for(plint j = domain.y0; j <= domain.y1; ++j)
			for(plint k = domain.z0; k <= domain.z1; ++k) {
				T * val = lattice.get(i, j, k).getExternal(Descriptor<T>::ExternalField::phaseBeginsAt);
				conc.get(i+offset.x, j+offset.y, k+offset.z) = *val;
			}
}

/*** ForceFieldFunctional ***/
template<class T, template<typename U> class Descriptor>
void ForceFieldFunctional<T, Descriptor>::process(
		Box3D domain,
		BlockLattice3D<T, Descriptor> & lattice,
		TensorField3D<T,3>& force)
{
	Dot3D offset = computeRelativeDisplacement(lattice, force);

	for(plint i = domain.x0; i <= domain.x1; ++i)
		for(plint j = domain.y0; j <= domain.y1; ++j)
			for(plint k = domain.z0; k <= domain.z1; ++k) {
				T * val = lattice.get(i, j, k).getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
				for(plint l=0; l < 3; ++l)
					force.get(i+offset.x, j+offset.y, k+offset.z)[l] = val[l];
			}
}

/*** ImmersedBoundaryVoxelizer3D ***/
template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryVoxelizer3D<T, Descriptor, Periodicity>::process(Box3D domain, BlockLattice3D<T, Descriptor> & lattice)
{
	if( ! ibm_.is_init()) {
		pcerr << "The immersed boundary module has not been initialized! Call the ImmersedBoundaryDynamics3D<...>::init() " <<
				"method after all particles have been added!" << std::endl;
		exit(-1);
	}

	/*Dot3D offset = lattice.getLocation();
	Box3D global_domain = domain.shift(offset.x, offset.y, offset.z);

	std::vector<Dot3D> points_inside;
	for(typename ImmersedBoundaryType::ObjMapIterator it = ibm_.particles_begin(); it != ibm_.particles_end(); ++it) {
		it->second->voxelizer().rehash();
		it->second->voxelizer().voxelize();
		it->second->voxelizer().to_dot_list(domain, points_inside);
	}

	// Set the phase of all points to 0
	for(plint i = domain.x0; i <= domain.x1; ++i)
		for(plint j = domain.y0; j <= domain.y1; ++j)
			for(plint k = domain.z0; k <= domain.z1; ++k)
				*(lattice.get(i, j, k).getExternal(Descriptor<T>::ExternalField::phaseBeginsAt)) = 0;

	// Insert
	for(plint i = 0; i < points_inside.size(); ++i) {
		T * phase = lattice.get(points_inside[i].x-offset.x, points_inside[i].y-offset.y, points_inside[i].z-offset.z)
				.getExternal(Descriptor<T>::ExternalField::phaseBeginsAt);
		*phase = 1.;
	}*/
}

template<class T, template<typename U> class Descriptor, class Periodicity>
ImmersedBoundaryVoxelizer3D<T, Descriptor, Periodicity> *
	wrap_ibm_voxelizer3D(ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity> & ibm)
{
	return new ImmersedBoundaryVoxelizer3D<T, Descriptor, Periodicity>(ibm);
}

/*** ImmersedBoundaryWrapperFunctional3D ***/
template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryWrapperFunctional3D<T, Descriptor, Periodicity>::process(
		Box3D domain,
		BlockLattice3D<T,Descriptor>& lattice,
		TensorField3D<T,3>& velocity)
{
	if( ! ibm.is_init()) {
		pcerr << "The immersed boundary module has not been initialized! Call the ImmersedBoundaryDynamics3D<...>::init() " <<
				"method after all particles have been added!" << std::endl;
		exit(-1);
	}

	ibm.interpolate_velocity(domain, velocity);
	ibm.compute_and_spread_forces(domain, lattice);
	ibm.move_vertices_and_revoxelize(domain, lattice);
}

/*template<class T, class Periodicity>
void ImmersedBoundaryWrapperFunctional3D<T, Descriptor, Periodicity>::process(
		Box3D domain,
		BlockLattice3D<T,Descriptor>& lattice,
		TensorField3D<T,3>& velocity)
{
	if( ! ibm.is_init()) {
		pcerr << "The immersed boundary module has not been initialized! Call the ImmersedBoundaryDynamics3D<...>::init() " <<
				"method after all particles have been added!" << std::endl;
		exit(-1);
	}

	ibm.interpolate_velocity(domain, velocity);
	ibm.compute_and_spread_forces(domain, lattice);
	ibm.move_vertices_and_revoxelize(domain, lattice);
}*/

template<class T, template<typename U> class Descriptor, class Periodicity>
ImmersedBoundaryWrapperFunctional3D<T, Descriptor, Periodicity> *
	wrap_ibm_dynamics3D(ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity> & ibm)
{
	return new ImmersedBoundaryWrapperFunctional3D<T, Descriptor, Periodicity>(ibm);
}

}

}



#endif /* DATAPROCESSORS_HH_ */
