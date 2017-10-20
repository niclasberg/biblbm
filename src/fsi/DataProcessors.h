/*
 * DataProcessors.h
 *
 *  Created on: Jun 24, 2015
 *      Author: niber
 */

#ifndef DATAPROCESSORS_H_
#define DATAPROCESSORS_H_
#include "ImmersedBoundaryDynamics.h"

namespace plb {

namespace fsi {

/*** ComputeConcentrationFunctional ***/
template<class T, template<typename U> class Descriptor>
class ComputeConcentrationFunctional : public BoxProcessingFunctional3D_LS<T, Descriptor, T> {
public:

	virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &, ScalarField3D<T> &);

	virtual ComputeConcentrationFunctional * clone() const
	{
		return new ComputeConcentrationFunctional(*this);
	}

	virtual BlockDomain::DomainT appliesTo() const {
		return BlockDomain::bulk;
	}

	virtual void getModificationPattern(std::vector<bool>& isWritten) const {
		isWritten[0] = false;
		isWritten[1] = true;
	}

	virtual void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
		modified[0] = modif::nothing;
		modified[1] = modif::staticVariables;
	}
};

/*** ForceFieldFunctional ***/
template<class T, template<typename U> class Descriptor>
class ForceFieldFunctional : public BoxProcessingFunctional3D_LT<T, Descriptor, T, 3> {
public:

	virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &, TensorField3D<T,3>&);

	virtual ForceFieldFunctional * clone() const
	{
		return new ForceFieldFunctional(*this);
	}

	virtual BlockDomain::DomainT appliesTo() const {
		return BlockDomain::bulk;
	}

	virtual void getModificationPattern(std::vector<bool>& isWritten) const {
		isWritten[0] = false;
		isWritten[1] = true;
	}

	virtual void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
		modified[0] = modif::nothing;
		modified[1] = modif::staticVariables;
	}
};

/*** ImmersedBoundaryVoxelizer3D ***/
template<class T, template<typename U> class Descriptor, class Periodicity>
class ImmersedBoundaryVoxelizer3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
	typedef ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity> ImmersedBoundaryType;

	ImmersedBoundaryVoxelizer3D(ImmersedBoundaryType & ibm)
	: ibm_(ibm)
	{ }

	virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> & field);

	virtual ImmersedBoundaryVoxelizer3D * clone() const
	{
		return new ImmersedBoundaryVoxelizer3D(*this);
	}

	virtual BlockDomain::DomainT appliesTo() const {
		return BlockDomain::bulk;
	}

	virtual void getModificationPattern(std::vector<bool>& isWritten) const {
		isWritten[0] = false;
	}

	virtual void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
		modified[0] = modif::staticVariables;
	}

private:
	ImmersedBoundaryType & ibm_;
};

template<class T, template<typename U> class Descriptor, class Periodicity>
ImmersedBoundaryVoxelizer3D<T, Descriptor, Periodicity> *
	wrap_ibm_voxelizer3D(ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity> & ibm);

// Wrapping of a general descriptor, no voxelization is done since we cannot know if 
// there is a phase field defined
template<class T, template<typename U> class Descriptor, class Periodicity>
class ImmersedBoundaryWrapperFunctional3D : public BoxProcessingFunctional3D_LT<T, Descriptor, T, 3> {
public:
	ImmersedBoundaryWrapperFunctional3D(
			ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity> & ibm_)
	: ibm(ibm_)
	{ }

	virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice, TensorField3D<T,3>& field);

	virtual ImmersedBoundaryWrapperFunctional3D * clone() const
	{
		return new ImmersedBoundaryWrapperFunctional3D(*this);
	}

	virtual BlockDomain::DomainT appliesTo() const {
		return BlockDomain::bulk;
	}

	virtual void getModificationPattern(std::vector<bool>& isWritten) const {
		isWritten[0] = true;
		isWritten[1] = false;
	}

	virtual void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
		modified[0] = modif::staticVariables;
		modified[1] = modif::nothing;
	}
private:
	ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity> & ibm;
};

// Specialization
/*template<class T, class Periodicity>
class ImmersedBoundaryWrapperFunctional3D : public BoxProcessingFunctional3D_LT<T, ForcedPhaseD3Q19Descriptor, T, 3> {
public:
	ImmersedBoundaryWrapperFunctional3D(
			ImmersedBoundaryDynamics3D<T, ForcedPhaseD3Q19Descriptor, Periodicity> & ibm_)
	: ibm(ibm_)
	{ }

	virtual void process(Box3D domain, BlockLattice3D<T,ForcedPhaseD3Q19Descriptor>& lattice, TensorField3D<T,3>& field);

	virtual ImmersedBoundaryWrapperFunctional3D * clone() const
	{
		return new ImmersedBoundaryWrapperFunctional3D(*this);
	}

	virtual BlockDomain::DomainT appliesTo() const {
		return BlockDomain::bulk;
	}

	virtual void getModificationPattern(std::vector<bool>& isWritten) const {
		isWritten[0] = true;
		isWritten[1] = false;
	}

	virtual void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
		modified[0] = modif::staticVariables;
		modified[1] = modif::nothing;
	}
private:
	ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity> & ibm;
};*/

template<class T, template<typename U> class Descriptor, class Periodicity>
ImmersedBoundaryWrapperFunctional3D<T, Descriptor, Periodicity> *
wrap_ibm_dynamics3D(ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity> & ibm);

}

}



#endif /* DATAPROCESSORS_H_ */
