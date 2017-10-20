#ifndef MACROPROCESSORS_HH_
#define MACROPROCESSORS_HH_
#include "MacroProcessors.h"

namespace plb {

template<class T, template<typename U> class Descriptor>
void update_density_velocity_force(
		MultiBlockLattice3D<T, Descriptor> & lattice,
		MultiScalarField3D<T> & density,
		MultiTensorField3D<T, Descriptor<T>::d> & velocity,
		Box3D domain,
		const Array<T, 3> & force)
{
	std::vector<MultiBlock3D *> blocks;
	blocks.push_back(&lattice);
	blocks.push_back(&density);
	blocks.push_back(&velocity);
	applyProcessingFunctional(new UpdateDensityVelocityFunctional<T, Descriptor>(force), domain, blocks);
}

template<class T, template<typename U> class Descriptor>
void UpdateDensityVelocityFunctional<T, Descriptor>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
	// Cast the blocks to the appropriate types
	BlockLattice3D<T, Descriptor> * lattice = dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[0]);
	ScalarField3D<T> * density = dynamic_cast<ScalarField3D<T> *>(blocks[1]);
	TensorField3D<T, Descriptor<T>::d> * velocity = dynamic_cast<TensorField3D<T, Descriptor<T>::d> *>(blocks[2]);

	Dot3D density_offset = computeRelativeDisplacement(*lattice, *density);
	Dot3D velocity_offset = computeRelativeDisplacement(*lattice, *velocity);
	Dot3D relative_position = lattice->getLocation();

	for(plint i = domain.x0; i <= domain.x1; ++i) {
		for(plint j = domain.y0; j <= domain.y1; ++j) {
			for(plint k = domain.z0; k <= domain.z1; ++k) {
				Cell<T, Descriptor> & cell = lattice->get(i, j, k);

				// Set the external force
				if(!boundary || boundary->contains(Array<T, 3>(relative_position.x + i, relative_position.y + j, relative_position.z + k))) {
					T * g = cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
					g[0] = external_force[0];
					g[1] = external_force[1];
					g[2] = external_force[2];
				}

				// Compute macroscopic variables
				density->get(i+density_offset.x, j+density_offset.y, k+density_offset.z) = cell.computeDensity();
				cell.computeVelocity(velocity->get(i+velocity_offset.x, j+velocity_offset.y, k+velocity_offset.z));
			}
		}
	}
}

template<class T, template<typename U> class Descriptor>
void update_velocity(
		MultiBlockLattice3D<T, Descriptor> & lattice,
		MultiTensorField3D<T, Descriptor<T>::d> & velocity,
		Box3D domain)
{
	std::vector<MultiBlock3D *> blocks;
	blocks.push_back(&lattice);
	blocks.push_back(&velocity);
	applyProcessingFunctional(new UpdateVelocityFunctional<T, Descriptor>, domain, blocks);
}

template<class T, template<typename U> class Descriptor>
void UpdateVelocityFunctional<T, Descriptor>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
	BlockLattice3D<T, Descriptor> * lattice = dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[0]);
	TensorField3D<T, Descriptor<T>::d> * velocity = dynamic_cast<TensorField3D<T, Descriptor<T>::d> *>(blocks[1]);
	Dot3D velocity_offset = computeRelativeDisplacement(*lattice, *velocity);

	for(plint i = domain.x0; i <= domain.x1; ++i) {
		for(plint j = domain.y0; j <= domain.y1; ++j) {
			for(plint k = domain.z0; k <= domain.z1; ++k) {
				// Compute total velocity
				lattice->get(i, j, k).computeVelocity(velocity->get(i+velocity_offset.x, j+velocity_offset.y, k+velocity_offset.z));
			}
		}
	}
}

}




#endif /* MACROPROCESSORS_HH_ */
