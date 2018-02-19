/*
 * deformation_test.cpp
 *
 *  Created on: Jun 2, 2015
 *      Author: niber
 */
#define FSI_PROFILE

#include "palabos3D.h"
#include "palabos3D.hh"
#include "fsi/headers.h"
#include "fsi/headers.hh"
#include "../common/calibrate_rbc_shape.h"
#include "../common/rbc_parameters.h"
#include "../common/shear_flow_setup.h"
//#include <random>

using namespace plb;
using namespace fsi;

#define DESCRIPTOR descriptors::ForcedPhaseD3Q19Descriptor
typedef Periodicity3D<double, true, false, true> Periodicity;

int main(int argc, char ** argv)
{
	plbInit(&argc, &argv);
	global::IOpolicy().activateParallelIO(true);
	global::directories().setOutputDir("test_dirac2/");
	global::directories().setInputDir(global::directories().getOutputDir());
	io::mkdir(global::directories().getOutputDir().c_str());

	plint margin_y = 1;
	plint geoNx = 36, geoNy = 38, geoNz = 36;
	Array<double, 3> channel_mid(0.5*geoNx, 0.5*geoNy + margin_y, 0.5*geoNz);
	plint Nt = 100000;
	plint output_interval = 100;

	// Number of particles in each direction
	plint np_x = 7;
	plint np_y = 2;
	plint np_z = 2;

	// Input parameters
	double tau = 1; 			// Relaxation parameter
	double Re_p = 0.1;			// Particle Reynolds number
	double lambda = 5;			// Viscosity ratio between the fluid inside and outside the rbcs
	double Ca = 0.1;			// Capillary number (mu*(rbc radius)*gamma / (rbc shear modulus)
	double k_bend_rel = 0.006;	// Relative bending energy
	double poisson_ratio = 0.990;// Poisson ratio

	// Load mesh
	ParticleShapeLibrary<double> shape_library;
	shape_library.read_and_store_mesh("sphere_TRImesh_626n_1248el.msh", "RBC");
	ParticleShape<double> * shape = shape_library.get_by_tag("RBC");
	double a_lb = shape->get_radius();

	// Deduce parameters
	double nu_lb = (tau - 0.5) / 3.0;
	double gamma_lb = Re_p *nu_lb / (a_lb*a_lb);
	double uMax = gamma_lb * geoNx / 2.0;

	// Relaxation frequencies
	double tau_inner = 0.5 + lambda*(tau - 0.5);
	double omega = 1. / tau;
	double omega_inner = 1. / tau_inner;

	// Membrane parameters
	double G_lb = a_lb * gamma_lb * nu_lb / Ca;
	double K_area = (poisson_ratio + 1)/(1 - poisson_ratio) * G_lb;
	double k_bend = k_bend_rel * G_lb * shape->get_area();
	RBCParameters<double> rbc_params = create_rbc_params(shape, G_lb,K_area, k_bend);

	RBCParticle<double> rbc(shape, rbc_params);

	pcout << "RBC Parameters:" << std::endl;
	pcout << "  Global area stiffness: " << rbc_params.k_area_global << std::endl;
	pcout << "  Local area stiffness: " << rbc_params.k_area_local << std::endl;
	pcout << "  Global volume stiffness: " << rbc_params.k_volume << std::endl;
	pcout << "  Bending stiffness: " << rbc_params.k_bend << std::endl;
	pcout << "  Equilibrium angle: " << rbc_params.theta0 << std::endl;
	pcout << "  Shear modulus: " << rbc_params.G() << std::endl;
	pcout << "  Area modulus " << rbc_params.K() << std::endl;
	pcout << "  Young's modulus: " << rbc_params.youngs_modulus() << std::endl;
	pcout << "==============" << std::endl;
	pcout << "Non-dimensional parameters:" << std::endl;
	pcout << "  Re_p: " << gamma_lb * util::sqr(shape->get_radius()) / nu_lb << std::endl;
	pcout << "  Ca_G: " << gamma_lb * shape->get_radius() * nu_lb / rbc_params.G() << std::endl;
	pcout << "  lambda: " << lambda << std::endl;
	pcout << "  Relative bending stiffness (kb/GA): " << rbc_params.k_bend / (rbc_params.G() * rbc.area()) << std::endl;
	pcout << "  Poisson ratio " << rbc_params.poisson_ratio() << std::endl;
	pcout << "==============" << std::endl;
	pcout << "Numerics:" << std::endl;
	pcout << "  Domain: " << geoNx << " x " << geoNy << " x " << geoNz << std::endl;
	pcout << "  uMax:" << uMax << std::endl;
	pcout << "  tau (outer fluid):" << tau << std::endl;
	pcout << "  tau (inner fluid):" << tau_inner << std::endl;
	pcout << "==============" << std::endl;

	// Create equilbrium shape
	shrink_rbc_volume(rbc, (double) 0.59*rbc.shape()->get_volume(), 10000);

	// Create lattice
	MultiBlockLattice3D<double, DESCRIPTOR> lattice(
			geoNx+1, geoNy+1 + 2*margin_y, geoNz+1,
			new GuoExternalForceVOFDynamics<double, DESCRIPTOR>(omega, omega_inner) );
	lattice.periodicity().toggle(0, true);
	lattice.periodicity().toggle(1, true);
	lattice.periodicity().toggle(2, true);
	Box3D domain = lattice.getBoundingBox();
	lattice.toggleInternalStatistics(false);

	OnLatticeBoundaryCondition3D<double,DESCRIPTOR> * boundaryCondition
			= createInterpBoundaryCondition3D<double,DESCRIPTOR>();
	shearFlowSetup(lattice, uMax, *boundaryCondition, margin_y);

	// Immersed boundary module
	ImmersedBoundaryDynamics3D<double, DESCRIPTOR, Periodicity> fsi(
			lattice,
			shape_library );

	// Domain boundary
	ParallelPlatesBoundary<double> boundary(margin_y, geoNy+margin_y);
	fsi.set_boundary(&boundary);
	//boundary.writeVTK("boundary.vtu");

	WallInteraction<double, MorsePotential<double> > wall_interaction(boundary, MorsePotential<double>(1, 1e-2, 1e-4));
	//fsi.add_force_decorator(&wall_interaction);

	// Velocity and concentration fields
	std::auto_ptr<MultiTensorField3D<double, 3> > velocity = generateMultiTensorField<double, 3>(domain);
	std::auto_ptr<MultiScalarField3D<double> > H = generateMultiScalarField<double>(domain);

	// Create rbcs
	double delta_x = (double) geoNx / (double)(np_x);
	double delta_y = (double) geoNy / (double)(np_y);
	double delta_z = (double) geoNz / (double)(np_z);

	rbc.set_minor_axis_orientation(Array<double, 3>(1, 0, 0));

	for(plint i = 0; i < np_x; ++i) {
		double x = ((double)i+0.5) * delta_x;
		for(plint j = 0; j< np_y; ++j) {
			double y = margin_y + ((double)j+0.5) * delta_y;
			for(plint k = 0; k < np_z; ++k) {
				double z = ((double)(k+0.5)) * delta_z;
				double randz[2] = {(double) std::rand() / (double) RAND_MAX - 0.5,
								   (double) std::rand() / (double) RAND_MAX - 0.5};
				global::mpi().bCast(randz, 2);

				Array<double, 3> pos = Array<double, 3>(x, y, z) + Array<double, 3>(0, randz[0] - 0.5*(y > (channel_mid[1]) ? -1.0 : 1.0), randz[1]);
				rbc.set_center_of_mass(pos);
				rbc.update();
				fsi.add_particle(&rbc);
			}
		}
	}

	// Initialize fsi
	fsi.init();

	pcout << "===============" << std::endl;
	pcout << "Starting coupled ib-lbm iterations" << std::endl;
	pcout << "Volume fraction: " << fsi.compute_total_volume() / (geoNx*geoNy*geoNz) << std::endl;

	for(plint it = 0; it <= Nt; ++it) {
		if((it % 1000) == 0) {
			pcout << "Iteration " << it << " / " << Nt << ", total cell volume: " << fsi.compute_total_volume() << std::endl;
		}

		computeVelocity(lattice, *velocity, domain);

		if((it % output_interval) == 0) {
			pcout << "Writing flow and particle data at iteration " << it << std::endl;
			fsi.write_particles_as_vtk(it);

			// Get concentration field
			applyProcessingFunctional(new ComputeConcentrationFunctional<double, DESCRIPTOR>, domain, lattice, *H);

			VtkImageOutput3D<double> vtkOut(createFileName("vtk", it, 6), 1.);
			vtkOut.writeData<3,float>(*velocity, "velocity", 1.);
			vtkOut.writeData<float>(*H, "concentration", 1.);
		}

		applyProcessingFunctional(new CouetteVelocityExtrapolator<double>(geoNy, margin_y, uMax), domain, *velocity);
		setExternalVector(lattice, domain, DESCRIPTOR<double>::ExternalField::forceBeginsAt, Array<double, 3>(0,0,0));

		applyProcessingFunctional(wrap_ibm_dynamics3D(fsi), domain, lattice, *velocity);

		// Compute forces and update node positions
		Profile::start_timer("lbm");
		lattice.collideAndStream();
		Profile::stop_timer("lbm");
	}

	Profile::write_report("profile");

	parallelIO::save(lattice, "checkpoint_lbm", false);
	fsi.save_checkpoint("checkpoint_fsi");
}

