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
#include "fsi/Time.h"
#include "../common/calibrate_rbc_shape.h"
#include "../common/channel_flow_setup.h"
//#include <random>

using namespace plb;
using namespace fsi;

#define DESCRIPTOR descriptors::ForcedPhaseD3Q19Descriptor
typedef Periodicity3D<double, true, false, true> Periodicity;

int main(int argc, char ** argv)
{
	plbInit(&argc, &argv);
	global::IOpolicy().activateParallelIO(true);
	global::directories().setOutputDir("Ca0_1_lambda_5_Rep0_1/");
	global::directories().setInputDir(global::directories().getOutputDir());
	io::mkdir(global::directories().getOutputDir().c_str());

	bool load_checkpoint = true;
	plint margin_y = 1;
	plint geoNx = 120, geoNy = 80, geoNz = 80;
	Array<double, 3> channel_mid(0.5*geoNx, 0.5*geoNy + margin_y, 0.5*geoNz);
	plint Nt = 600000;
	plint output_interval = 500;

	// Number of particles in each direction
	plint np_x = 20;
	plint np_y = 5;
	plint np_z = 4;

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
	double uMax = gamma_lb * geoNy / 4.0;

	// Relaxation frequencies
	double tau_inner = 0.5 + lambda*(tau - 0.5);
	double omega = 1. / tau;
	double omega_inner = 1. / tau_inner;

	// Body force
	Array<double, 3> force(2*gamma_lb*nu_lb / ((double) geoNy), 0, 0);

	// Membrane parameters
	RBCParameters<double> rbc_params;
	double G_lb = a_lb * gamma_lb * nu_lb / Ca;

	// In plane properties
	// Get average link length
	double lavg = 0;
	for(ParticleShape<double>::link_const_iterator it = shape->links_begin(); it != shape->links_end(); ++it)
		lavg += it->length;
	lavg /= shape->count_links();

	rbc_params.p = 1;
	rbc_params.L0 = lavg;
	rbc_params.vol_desired = shape->get_volume();
	rbc_params.Lmax = 2.2*rbc_params.L0;
	rbc_params.k_in_plane = 1;
	rbc_params.k_in_plane = G_lb / rbc_params.G();
	rbc_params.compute_C();

	// Bending stiffness
	rbc_params.k_bend = k_bend_rel * G_lb * shape->get_area();
	rbc_params.theta0 = std::acos((std::sqrt(3)*(shape->count_vertices() - 2) - 5*M_PI) / (std::sqrt(3)*(shape->count_vertices() - 2) - 3*M_PI) );
	//rbc_params.theta0 = 0;

	// Area and volume constraint
	double K_area = (poisson_ratio + 1)/(1 - poisson_ratio) * G_lb;
	double delta_K_area = K_area - rbc_params.K();

	rbc_params.k_area_global = delta_K_area * 0.9;
	rbc_params.k_area_local = delta_K_area * 0.1;
	rbc_params.k_volume = 0.5*delta_K_area;

	RBCParticle<double> rbc(shape_library.get_by_tag("RBC"), rbc_params);

	pcout << "RBC Parameters:" << std::endl;
	pcout << "  Global area stiffness: " << rbc_params.k_area_global << std::endl;
	pcout << "  Local area stiffness: " << rbc_params.k_area_local << std::endl;
	pcout << "  Global volume stiffness: " << rbc_params.k_volume << std::endl;
	pcout << "  Bending stiffness: " << rbc_params.k_bend << std::endl;
	pcout << "  In-plane shear stiffness: " << rbc_params.k_in_plane << std::endl;
	pcout << "  Equilibrium angle: " << rbc_params.theta0 << std::endl;
	pcout << "  Shear modulus: " << rbc_params.G() << std::endl;
	pcout << "  Area modulus " << rbc_params.K() << std::endl;
	pcout << "  Young's modulus: " << rbc_params.youngs_modulus() << std::endl;
	pcout << "  Area dilation energy " << rbc_params.C << std::endl;
	pcout << "==============" << std::endl;
	pcout << "Non-dimensional parameters:" << std::endl;
	pcout << "  Re_H: " << (uMax/2.0) * geoNy / nu_lb << std::endl;
	pcout << "  Re_p: " << gamma_lb * util::sqr(shape->get_radius()) / nu_lb << std::endl;
	pcout << "  Ca_G: " << gamma_lb * shape->get_radius() * nu_lb / rbc_params.G() << std::endl;
	pcout << "  lambda: " << lambda << std::endl;
	pcout << "  Relative bending stiffness (kb/GA): " << rbc_params.k_bend / (rbc_params.G() * rbc.area()) << std::endl;
	pcout << "  Poisson ratio: " << rbc_params.poisson_ratio() << std::endl;
	pcout << "==============" << std::endl;
	pcout << "Numerics:" << std::endl;
	pcout << "  Domain: " << geoNx << " x " << geoNy << " x " << geoNz << std::endl;
	pcout << "  uMax:" << uMax << std::endl;
	pcout << "  tau (outer fluid):" << tau << std::endl;
	pcout << "  tau (inner fluid):" << tau_inner << std::endl;
	pcout << "==============" << std::endl;

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
	channelFlowSetup(lattice, uMax, *boundaryCondition);

	if(load_checkpoint)
		parallelIO::load("checkpoint_lbm", lattice, false);

	// Immersed boundary module
	ImmersedBoundaryDynamics3D<double, DESCRIPTOR, Periodicity> fsi(
			lattice,
			shape_library );

	// Domain boundary
	ParallelPlatesBoundary<double> boundary(margin_y, geoNy+margin_y);
	fsi.set_boundary(&boundary);

	// Create rbcs
	if(load_checkpoint) {
		fsi.load_checkpoint("checkpoint_fsi");
	} else {
		// Create equilbrium shape
		shrink_rbc_volume(rbc, (double) 0.59*rbc.shape()->get_volume(), 10000);

		double y0 = 2. + margin_y;
		double delta_x = (double) geoNx / (double)(np_x);
		double delta_y = (double) (geoNy-4) / (double)(np_y);
		double delta_z = (double) geoNz / (double)(np_z);

		rbc.set_minor_axis_orientation(Array<double, 3>(1, 0, 0));

		for(plint i = 0; i < np_x; ++i) {
			double x = ((double)i+0.5) * delta_x;
			for(plint j = 0; j< np_y; ++j) {
				double y = y0 + ((double)j+0.5) * delta_y;
				for(plint k = 0; k < np_z; ++k) {
					double z = ((double)(k+0.5) + ((j%2==0) ? 0 : 0.5)) * delta_z;
					double randz[3] = {(double) std::rand() / (double) RAND_MAX - 0.5,
									   (double) std::rand() / (double) RAND_MAX - 0.5,
									   (double) std::rand() / (double) RAND_MAX - 0.5};
					global::mpi().bCast(randz, 3);

					Array<double, 3> pos = Array<double, 3>(x, y, z) + Array<double, 3>(1.0*randz[0], 0.5*randz[1], randz[2]);
					rbc.set_center_of_mass(pos);
					rbc.update();
					fsi.add_particle(&rbc);
				}
			}
		}
	}

	// Set iteration
	Time it(0);
	if(load_checkpoint) {
		it.read_checkpoint("checkpoint_time");
	}

	// Velocity and concentration fields
	std::auto_ptr<MultiTensorField3D<double, 3> > velocity = generateMultiTensorField<double, 3>(domain);
	std::auto_ptr<MultiScalarField3D<double> > H = generateMultiScalarField<double>(domain);

	// Initialize fsi
	fsi.init();

	pcout << "===============" << std::endl;
	pcout << "Starting coupled ib-lbm iterations" << std::endl;
	pcout << "Volume fraction: " << fsi.compute_total_volume() / (geoNx*geoNy*geoNz) << std::endl;

	for( ; it <= Nt; ++it) {
		if((it.to_plint() % 1) == 0) {
			pcout << "Iteration " << it << " / " << Nt << ", total cell volume: " << fsi.compute_total_volume() << std::endl;
		}

		computeVelocity(lattice, *velocity, domain);
		//applyProcessingFunctional(new ChannelVelocityExtrapolator<double>(geoNy, margin_y, uMax), domain, *velocity);
		setExternalVector(lattice, domain, DESCRIPTOR<double>::ExternalField::forceBeginsAt, force);

		applyProcessingFunctional(wrap_ibm_dynamics3D(fsi), domain, lattice, *velocity);

		if((it.to_plint() % output_interval) == 0) {
			fsi.write_particles_as_vtk(it.to_plint());

			// Get concentration field
			applyProcessingFunctional(new ComputeConcentrationFunctional<double, DESCRIPTOR>, domain, lattice, *H);

			VtkImageOutput3D<double> vtkOut(createFileName("vtk", it.to_plint(), 6), 1.);
			vtkOut.writeData<3,float>(*velocity, "velocity", 1.);
			vtkOut.writeData<float>(*H, "concentration", 1.);
		}

		// Compute forces and update node positions
		Profile::start_timer("lbm");
		lattice.collideAndStream();
		Profile::stop_timer("lbm");
	}

	Profile::write_report("profile");

	parallelIO::save(lattice, "checkpoint_lbm", false);
	fsi.save_checkpoint("checkpoint_fsi");
	it.write_checkpoint("checkpoint_time");
}

