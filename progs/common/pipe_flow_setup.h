/*
 * pipe_flow_setup.h
 *
 *  Created on: Jun 30, 2015
 *      Author: niber
 */

#ifndef PIPE_FLOW_SETUP_H_
#define PIPE_FLOW_SETUP_H_

template<class T>
struct PipeFlowParameters {

};

template<class T>
void read_pipe_flow_parameters(std::string fname, plb::fsi::ParticleShapeLibrary<T> & shape_library, PipeFlowParameters<T> & params)
{

}

/* Flow initialization */
template<class T>
class PipeFlowDensityAndVelocity {
public:
	PipeFlowDensityAndVelocity(const plb::fsi::PipeBoundary<T> & boundary_, T lattice_u_)
	: boundary(boundary_), lattice_u(lattice_u_)
	{

	}

	void operator()(plb::plint iX, plb::plint iY, plb::plint iZ, T & rho, plb::Array<T, 3> & u) const
	{
		rho = (T) 1;
		T r_sqr = (plb::util::sqr(iY - boundary.y0) + plb::util::sqr(iZ - boundary.z0)) / plb::util::sqr(boundary.radius);
		u[0] = (T) 0; //lattice_u * ((T) 1 - r_sqr);
		u[1] = (T) 0;
		u[2] = (T) 0;
	}
private:
	const plb::fsi::PipeBoundary<T> & boundary;
	T lattice_u;
};



template<class T, template<typename U> class Descriptor>
void pipeFlowSetup(
		plb::MultiBlockLattice3D<T, Descriptor> & lattice,
		const plb::Array<T, 3> & u_max,
		const plb::Array<T, 3> & external_force,
		plb::fsi::PipeBoundary<T> * boundary)
{
	lattice.periodicity().toggle(0, true);
	lattice.periodicity().toggle(1, false);
	lattice.periodicity().toggle(2, false);

	// Setup boundary condition for LBM
	defineDynamics<double,Descriptor>(lattice,
			lattice.getBoundingBox(),
			new BoundarySurface<T>(*boundary),
			new plb::BounceBack<T,Descriptor>() );
	defineDynamics<double,Descriptor>(lattice,
			lattice.getBoundingBox(),
			new BoundaryOutside<T>(*boundary),
			new plb::NoDynamics<T, Descriptor>() );

	// Set external force
	setExternalVector(lattice, lattice.getBoundingBox(), Descriptor<T>::ExternalField::forceBeginsAt, external_force);

	// Initialize flow field
	initializeAtEquilibrium(lattice,
			lattice.getBoundingBox(),
			PipeFlowDensityAndVelocity<T>(*boundary, u_max) );
}



#endif /* PIPE_FLOW_SETUP_H_ */
