/*
 * MacroProcessors.h
 *
 *  Created on: Feb 25, 2014
 *      Author: niber
 */

#ifndef MACROPROCESSORS_H_
#define MACROPROCESSORS_H_
#include "Boundary.h"

namespace plb {

template<class T, template<typename U> class Descriptor>
void update_density_velocity_force(
		MultiBlockLattice3D<T, Descriptor> &,
		MultiScalarField3D<T> &,
		MultiTensorField3D<T, Descriptor<T>::d> &,
		Box3D,
		const Array<T, 3> &);

template<class T, template<typename U> class Descriptor>
class UpdateDensityVelocityFunctional : public BoxProcessingFunctional3D {
public:
	UpdateDensityVelocityFunctional(const Array<T, 3> & force) : external_force(force), boundary(0) { }
	UpdateDensityVelocityFunctional(const Array<T, 3> & force, fsi::Boundary<T> * _boundary)
	: external_force(force), boundary(_boundary) { }

	virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *>);

	virtual UpdateDensityVelocityFunctional * clone() const
	{
		return new UpdateDensityVelocityFunctional(*this);
	}

	virtual BlockDomain::DomainT appliesTo() const {
		return BlockDomain::bulk;
	}

	virtual void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
		// Note: The 2nd and 3rd fields (density and velocity) are actually written,
		// but we put the modification to "nothing" in order to suppress the communication.
		// Only the bulk field is needed for the fsi.
		modified[0] = modif::nothing;
		modified[1] = modif::nothing;
		modified[2] = modif::nothing;
	}

private:
	fsi::Boundary<T> * boundary;
	Array<T, 3> external_force;
};

template<class T, template<typename U> class Descriptor>
void update_velocity(
		MultiBlockLattice3D<T, Descriptor> & lattice,
		MultiTensorField3D<T, Descriptor<T>::d> & velocity,
		Box3D domain);

template<class T, template<typename U> class Descriptor>
class UpdateVelocityFunctional : public BoxProcessingFunctional3D {
public:
	virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *>);

	virtual UpdateVelocityFunctional * clone() const
	{
		return new UpdateVelocityFunctional(*this);
	}

	virtual BlockDomain::DomainT appliesTo() const {
		return BlockDomain::bulk;
	}

	virtual void getModificationPattern(std::vector<bool>& isWritten) const {
		isWritten[0] = false;
		isWritten[1] = false;
	}

	virtual void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
		modified[0] = modif::nothing;
		modified[1] = modif::nothing;
	}
};

}



#endif /* MACROPROCESSORS_H_ */
