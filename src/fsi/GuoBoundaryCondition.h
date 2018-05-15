/*
 * GuoBoundaryCondition.h
 *
 *  Created on: Aug 17, 2015
 *      Author: niber
 */

#ifndef GUOBOUNDARYCONDITION_H_
#define GUOBOUNDARYCONDITION_H_
#include <vector>

namespace plb {

namespace fsi {

template<class T> class Boundary;

template<class T, template<typename U> class Descriptor>
class GuoRigidWallBoundary {
public:
	struct NodeInfo {
		Dot3D location;
		T delta;
		Dot3D direction;
		Array<T, 3> wall_position;
		Array<T, 3> velocity;
	};

	typedef typename std::vector<NodeInfo>::iterator iterator;
	typedef typename std::vector<NodeInfo>::const_iterator const_iterator;

	GuoRigidWallBoundary(const Boundary<T> & boundary, MultiBlockLattice3D<T, Descriptor> & lattice);
	const Boundary<T> & boundary() const { return boundary_; }

	void insert();
	void put_node(Dot3D pos, Dot3D dir, T dist, const Array<T, 3> &);

	void set_wall_velocity(const Array<T, 3> & velocity);
	template<class VelocityFunction> void set_wall_velocity(VelocityFunction velocity_function);

	iterator begin() { return wall_nodes_.begin(); }
	const_iterator begin() const { return wall_nodes_.begin(); }
	iterator end() { return wall_nodes_.end(); }
	const_iterator end() const { return wall_nodes_.end(); }


private:
	std::vector<NodeInfo> wall_nodes_;
	const Boundary<T> & boundary_;
	MultiBlockLattice3D<T, Descriptor> & lattice_;
};

template<class T, template<typename U> class Descriptor>
class GuoRigidWallBoundaryInstantiator : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
	GuoRigidWallBoundaryInstantiator(GuoRigidWallBoundary<T, Descriptor> & boundary_model)
	: model_(boundary_model)
	{

	}

	virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &);

	virtual GuoRigidWallBoundaryInstantiator * clone() const
	{
		return new GuoRigidWallBoundaryInstantiator(*this);
	}

	virtual BlockDomain::DomainT appliesTo() const {
		return BlockDomain::bulk;
	}

	virtual void getModificationPattern(std::vector<bool>& isWritten) const {
		isWritten[0] = false;
	}

	virtual void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
		modified[0] = modif::nothing;
	}
private:
	GuoRigidWallBoundary<T, Descriptor> & model_;
};


template<class T, template<typename U> class Descriptor>
class GuoRigidWallBoundaryFunctional : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
	GuoRigidWallBoundaryFunctional(GuoRigidWallBoundary<T, Descriptor> & model) : model_(model) { }

	virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &);

	virtual GuoRigidWallBoundaryFunctional * clone() const
	{
		return new GuoRigidWallBoundaryFunctional(*this);
	}

	virtual BlockDomain::DomainT appliesTo() const {
		return BlockDomain::bulk;
	}

	virtual void getModificationPattern(std::vector<bool>& isWritten) const {
		isWritten[0] = true;
	}

	virtual void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
		modified[0] = modif::staticVariables;
	}

private:
	GuoRigidWallBoundary<T, Descriptor> & model_;
};

template<class T, template<typename U> class Descriptor>
class GuoRigidWallBoundaryDebugger : public BoxProcessingFunctional3D_S<T> {
public:
	GuoRigidWallBoundaryDebugger(GuoRigidWallBoundary<T, Descriptor> & model) : model_(model) { }

	virtual void process(Box3D domain, ScalarField3D<T> &);

	virtual GuoRigidWallBoundaryDebugger * clone() const
	{
		return new GuoRigidWallBoundaryDebugger(*this);
	}

	virtual BlockDomain::DomainT appliesTo() const {
		return BlockDomain::bulk;
	}

	virtual void getModificationPattern(std::vector<bool>& isWritten) const {
		isWritten[0] = true;
	}

	virtual void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
		modified[0] = modif::staticVariables;
	}

private:
	GuoRigidWallBoundary<T, Descriptor> & model_;
};

}

}



#endif /* GUOBOUNDARYCONDITION_H_ */
