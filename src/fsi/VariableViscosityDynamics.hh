/*
 * VariableViscosityDynamics.hh
 *
 *  Created on: Jun 24, 2015
 *      Author: niber
 */

#ifndef VARIABLEVISCOSITYDYNAMICS_HH_
#define VARIABLEVISCOSITYDYNAMICS_HH_
#include "VariableViscosityDynamics.h"

namespace plb {

/*** ForcedPhaseD3Q19Descriptor ***/
namespace descriptors {

template<class T>
const char ForcedPhaseD3Q19Descriptor<T>::name[] = "ForcedPhaseD3Q19";

} /* namespace descriptors */


/*
 * Now, since we've defined a new descriptor, the optimized force evaluation that
 * is defined for the ForcedD3Q19Descriptor will no longer be called. We hence
 * redefine these methods. We also introduce a helper function to compute the
 * omega (linear interpolation of viscosity)
 */
template<typename T, typename Descriptor>
struct vofTemplatesImpl
{
	static void addGuoForce(
			Array<T,Descriptor::q>& f,
			T* externalScalars,
			Array<T,Descriptor::d> const& u, T omega, T amplitude );

	static void computeOmega(
			T * externalScalars,
			T omega1,
			T omega2,
			T & omega
			);
};

template<typename T>
struct vofTemplatesImpl<T, descriptors::ForcedPhaseD3Q19Descriptor<T> > {
	static void addGuoForce(
					Array<T,descriptors::ForcedPhaseD3Q19Descriptor<T>::q>& f,
					T* externalScalars,
					Array<T,descriptors::ForcedPhaseD3Q19Descriptor<T>::d> const& u, T omega, T amplitude )
	{
		static const int forceBeginsAt
			= descriptors::ForcedPhaseD3Q19Descriptor<T>::ExternalField::forceBeginsAt;
		T* force = externalScalars + forceBeginsAt;
		T mu = amplitude*((T)1-omega/(T)2);

		static const T oneOver6 = (T)1/(T)6;
		static const T oneOver12 = (T)1/(T)12;

		f[0] += -mu*(force[0]*u[0]+force[1]*u[1]+force[2]*u[2]);

		f[1] += oneOver6*mu*(force[0]*(-(T)1+2*u[0])-force[1]*u[1]-force[2]*u[2]);

		f[2] += -oneOver6*mu*(force[0]*u[0]+force[1]*((T)1-2*u[1])+force[2]*u[2]);

		f[3] += -oneOver6*mu*(force[0]*u[0]
							   + force[1]*u[1]
							   + force[2]*((T)1-2*u[2]));

		f[4] += oneOver12*mu*(force[0]*(-(T)1+2*u[0]+3*u[1])
							   + force[1]*(-(T)1+2*u[1]+3*u[0])
							   - force[2]*u[2]);

		f[5] += oneOver12*mu*( force[0]*(-(T)1+2*u[0]-3*u[1])
								+ force[1]*((T)1+2*u[1]-3*u[0])
								- force[2]*u[2]);

		f[6] += oneOver12*mu*(force[0]*(-(T)1+2*u[0]+3*u[2])
							   - force[1]*u[1]
							   + force[2]*(-(T)1+2*u[2]+3*u[0]));

		f[7] += oneOver12*mu*(force[0]*(-(T)1+2*u[0]-3*u[2])
							   - force[1]*u[1]
							   + force[2]*((T)1+2*u[2]-3*u[0]));

		f[8] += -oneOver12*mu*(force[0]*u[0]
								+ force[1]*((T)1-2*u[1]-3*u[2])
								+ force[2]*((T)1-2*u[2]-3*u[1]));

		f[9] += -oneOver12*mu*(force[0]*u[0]
								+ force[1]*((T)1-2*u[1]+3*u[2])
								+ force[2]*(-(T)1-2*u[2]+3*u[1]));

		f[10] += oneOver6*mu*(force[0]*((T)1+2*u[0])
								-force[1]*u[1]
								-force[2]*u[2]);

		f[11] += -oneOver6*mu*(force[0]*u[0]
								+force[1]*(-(T)1-2*u[1])
								+force[2]*u[2]);

		f[12] += -oneOver6*mu*(force[0]*u[0]
								 +force[1]*u[1]
								 +force[2]*(-(T)1-2*u[2]));

		f[13] += oneOver12*mu*(force[0]*((T)1+2*u[0]+3*u[1])
								 +force[1]*((T)1+2*u[1]+3*u[0])
								 -force[2]*u[2]);

		f[14] += oneOver12*mu*(force[0]*((T)1+2*u[0]-3*u[1])
								 +force[1]*(-(T)1+2*u[1]-3*u[0])
								 -force[2]*u[2]);

		f[15] += oneOver12*mu*(force[0]*((T)1+2*u[0]+3*u[2])
								 -force[1]*u[1]
								 +force[2]*((T)1+2*u[2]+3*u[0]));

		f[16] += oneOver12*mu*(force[0]*((T)1+2*u[0]-3*u[2])
								 -force[1]*u[1]
								 +force[2]*(-(T)1+2*u[2]-3*u[0]));

		f[17] += -oneOver12*mu*(force[0]*u[0]
								  +force[1]*(-(T)1-2*u[1]-3*u[2])
								  +force[2]*(-(T)1-2*u[2]-3*u[1]));

		f[18] += -oneOver12*mu*(force[0]*u[0]
								  +force[1]*(-(T)1-2*u[1]+3*u[2])
								  +force[2]*((T)1-2*u[2]+3*u[1]));
	}

	static void computeOmega(T * externalScalars, T omega0, T omega1, T & omega)
	{
		static const int phaseBeginsAt
					= descriptors::ForcedPhaseD3Q19Descriptor<T>::ExternalField::phaseBeginsAt;
		const T c = externalScalars[phaseBeginsAt];

		// Linear dependence on concentration:
		//   nu = c*nu_0 + (1-c)*nu_1
		//   where nu_i = cs^2 * (1/omega_i - 1/2)
		omega = omega0*omega1 / (c*omega0 + ((T)1.-c)*omega1);
	}

};

/*** GuoExternalForceVOFDynamics ***/
template<typename T, template<typename U> class Descriptor>
int GuoExternalForceVOFDynamics<T, Descriptor>::id =
		meta::registerTwoParamDynamics<T, Descriptor, GuoExternalForceVOFDynamics<T, Descriptor> >("Guo_VOF_dynamics");

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceVOFDynamics<T, Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
	T omegaTot;
	T rhoBar = this->computeRhoBar(cell);
    Array<T,Descriptor<T>::d> u, j;
    this->computeVelocity(cell, u);
    T rho = Descriptor<T>::fullRho(rhoBar);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
        j[iD] = rho * u[iD];

    vofTemplatesImpl<T, Descriptor<T> >::computeOmega(cell.getExternal(0), this->getOmega(), this->omega2_, omegaTot);

    T uSqr = dynamicsTemplates<T,Descriptor>::bgk_ma2_collision(cell, rhoBar, j, omegaTot);
    //externalForceTemplates<T,Descriptor>::addGuoForce(cell, u, omegaTot, (T) 1);
    vofTemplatesImpl<T,Descriptor<T> >::addGuoForce(cell.getRawPopulations(), cell.getExternal(0), u, omegaTot, (T)1);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T GuoExternalForceVOFDynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}


template<typename T, template<typename U> class Descriptor>
GuoExternalForceVOFDynamics<T,Descriptor>* GuoExternalForceVOFDynamics<T,Descriptor>::clone() const
{
	return new GuoExternalForceVOFDynamics<T, Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int GuoExternalForceVOFDynamics<T,Descriptor>::getId() const
{
	return id;
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceVOFDynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
	serializer.addValue(this->getOmega());
	serializer.addValue(omega2_);
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceVOFDynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
	this->setOmega(unserializer.readValue<T>());
	unserializer.readValue(omega2_);
}


}


#endif /* VARIABLEVISCOSITYDYNAMICS_HH_ */
