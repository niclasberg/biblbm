#ifndef EXTENSIONAL_FLOW_SETUP_H_
#define EXTENSIONAL_FLOW_SETUP_H_

////// ________________________________________________________________
template<class T>
inline plb::Array<T, 3> extensionalFlowVelocity(plb::plint iX, plb::plint iY, plb::plint geoNx, plb::plint geoNy, T eps) {
	const T x = iX - (T)geoNx / 2.;
	const T y = iY - (T)geoNy / 2.;
	return plb::Array<T, 3>(eps*x, -eps*y, 0);
}

////// ________________________________________________________________
template <typename T>
class ExtensionalFlowVelocity {
public:
	ExtensionalFlowVelocity(plb::plint geoNx_, plb::plint geoNy_, T eps_)
	: geoNx(geoNx_), geoNy(geoNy_), eps(eps_)
	{ }
	void operator()(plb::plint iX, plb::plint iY, plb::plint iZ, plb::Array<T,3>& u) const  {
		u = extensionalFlowVelocity(iX, iY, geoNx, geoNy, eps);
	}
private:
	plb::plint geoNx;
	plb::plint geoNy;
	T eps;
};

////// ________________________________________________________________
template <typename T>
class ExtensionalFlowDensityAndVelocity {
public:
	ExtensionalFlowDensityAndVelocity(plb::plint geoNx_, plb::plint geoNy_, T eps_)
	: geoNx(geoNx_), geoNy(geoNy_), eps(eps_)
	{ }
	void operator()(plb::plint iX, plb::plint iY, plb::plint iZ, T &rho, plb::Array<T,3>& u) const {
		rho = (T)1;
		u = extensionalFlowVelocity(iX, iY, geoNx, geoNy, eps);
	}
private:
	plb::plint geoNx;
	plb::plint geoNy;
	T eps;
};

////// ________________________________________________________________
template<class T, template<typename U> class Descriptor>
void extensionalFlowSetup( 
		plb::MultiBlockLattice3D<T,Descriptor>& lattice,
		T eps,
		plb::OnLatticeBoundaryCondition3D<T,Descriptor>& boundaryCondition)
{
	using namespace plb;

	const plint nx = lattice.getBoundingBox().getNx();
	const plint ny = lattice.getBoundingBox().getNy();
	const plint nz = lattice.getBoundingBox().getNz();
	const plint geoNx = nx-1;
	const plint geoNy = ny-1;
	const plint geoNz = nz-1;

	plb::Box3D leftSideSurface (0,     0,       1,     geoNy-1, 0, geoNz);
	plb::Box3D rightSideSurface(geoNx, geoNx,   1,     geoNy-1, 0, geoNz);
	plb::Box3D topSurface      (1,     geoNx-1, geoNy, geoNy,   0, geoNz);
	plb::Box3D bottomSurface   (1,     geoNx-1, 0,     0,       0, geoNz);

	// Corner edges
	plb::Box3D bottomLeftEdge (0,     0,     0,     0,     0, geoNz);
	plb::Box3D bottomRightEdge(geoNx, geoNx, 0,     0,     0, geoNz);
	plb::Box3D topLeftEdge    (0,     0,     geoNy, geoNy, 0, geoNz);
	plb::Box3D topRightEdge   (geoNx, geoNx, geoNy, geoNy, 0, geoNz);

	lattice.periodicity().toggle(0, false); 
	lattice.periodicity().toggle(1, false); 
	lattice.periodicity().toggle(2, true); // Use periodic boundaries on Z direction.

	// Add boundary conditions along the left, right, top and bottom walls
	boundaryCondition.addVelocityBoundary0N(leftSideSurface, lattice);
	boundaryCondition.addVelocityBoundary0P(rightSideSurface,    lattice);
	boundaryCondition.addVelocityBoundary1N(bottomSurface, lattice);
	boundaryCondition.addVelocityBoundary1P(topSurface,    lattice);

	// Corner conditions
	boundaryCondition.addExternalVelocityEdge2NN(bottomLeftEdge, lattice);
	boundaryCondition.addExternalVelocityEdge2PN(bottomRightEdge, lattice);
	boundaryCondition.addExternalVelocityEdge2NP(topLeftEdge, lattice);
	boundaryCondition.addExternalVelocityEdge2PP(topRightEdge, lattice);

	initializeAtEquilibrium(lattice, lattice.getBoundingBox(),
			ExtensionalFlowDensityAndVelocity<T>(geoNx, geoNy, eps) );

	setBoundaryVelocity(lattice, lattice.getBoundingBox(), ExtensionalFlowVelocity<T>(geoNx, geoNy, eps) );

	lattice.initialize();
}


#endif /* EXTENSIONAL_FLOW_SETUP_H_ */
