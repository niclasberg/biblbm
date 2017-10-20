/*
 * VariableViscosityDynamics.h
 *
 *  Created on: Jun 23, 2015
 *      Author: niber
 */

#ifndef VARIABLEVISCOSITYDYNAMICS_H_
#define VARIABLEVISCOSITYDYNAMICS_H_

#include "core/globalDefs.h"
#include "core/dynamics.h"
#include "basicDynamics/isoThermalDynamics.h"

namespace plb {

/*
 * In whole blood, the viscosity of the fluid inside
 * the RBCs is different than that outside. The phase
 * information will be stored in the external field of
 * each cell. Hence a new descriptor needs to be defined
 * to hold this information
 */
namespace descriptors {
struct ForcedPhasedescriptor3D {
    static const int numScalars     = 4;
    static const int phaseBeginsAt = 3;
    static const int forceBeginsAt  = 0;
    static const int sizeOfForce    = 3;
};

struct ForcedPhasedescriptorBase3D {
    typedef ForcedPhasedescriptor3D ExternalField;
};

// To make life easier, we only implement a descriptor for the D3Q19 lattice
template <typename T>
struct ForcedPhaseD3Q19Descriptor
: public D3Q19DescriptorBase<T>, public ForcedPhasedescriptorBase3D
{
	static const char name[];
};

} /* namespace descriptors */

// Dynamics object implementing the variable viscosity
template<class T, template<typename U> class Descriptor>
class GuoExternalForceVOFDynamics : public  ExternalForceDynamics<T,Descriptor> {
public:
	GuoExternalForceVOFDynamics(T omega1, T omega2)
	:  ExternalForceDynamics<T,Descriptor>(omega1), omega2_(omega2)
	  { }

    /// Clone the object on its dynamic type.
    GuoExternalForceVOFDynamics<T,Descriptor>* clone() const;

    /// Implementation of the collision step
	virtual void collide(Cell<T,Descriptor>& cell,
						 BlockStatistics& statistics_);

	/// Compute equilibrium distribution function
	virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
								 T jSqr, T thetaBar=T()) const;

	/// Return a unique ID for this class.
	virtual int getId() const;

	/// Serialize the dynamics object.
	virtual void serialize(HierarchicSerializer& serializer) const;
	virtual void unserialize(HierarchicUnserializer& unserializer);

private:
    T omega2_;
    static int id;
};

}



#endif /* VARIABLEVISCOSITYDYNAMICS_H_ */
