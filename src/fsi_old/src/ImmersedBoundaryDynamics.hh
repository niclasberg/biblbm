#ifndef IMMERSEDBOUNDARYDYNAMICS_HH_
#define IMMERSEDBOUNDARYDYNAMICS_HH_
#include "ImmersedBoundaryDynamics.h"
#include "atomicBlock/dataProcessor3D.h"

namespace plb {

template<class T, template<typename U> class Descriptor>
class FsiForceProcessor : public BoxProcessingFunctional3D_
{
public:

private:

};

template<class T, template<typename U> class Descriptor>
void ImmersedBoundaryDynamics3D<T, Descriptor>::compute_fsi_forces(
		const MultiBlockLattice3D<T, Descriptor> & lattice,
		const MultiTensorField3D<T, 3> & velocity
)
{
	//
}

}




#endif /* IMMERSEDBOUNDARYDYNAMICS_HH_ */
