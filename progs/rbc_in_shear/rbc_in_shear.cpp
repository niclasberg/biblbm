// Comment away this line to disable profiling:
//#define FSI_PROFILE

#include "palabos3D.h"
#include "palabos3D.hh"
#include "fsi/headers.h"
#include "fsi/headers.hh"
#include "../common/calibrate_rbc_shape.h"
#include "../common/shear_flow_setup.h"
#include "../common/rbc_parameters.h"

using namespace plb;
using namespace fsi;

#define DESCRIPTOR descriptors::ForcedPhaseD3Q19Descriptor

// Define the periodicity of the domain (periodicity_x, periodicity_y, periodicity_z)
typedef Periodicity3D<double, true, false, true> Periodicity;

int main(int argc, char ** argv)
{
	// Initialize palabos
	plbInit(&argc, &argv);

	// Setup input/output
	global::IOpolicy().activateParallelIO(true);

	// Domain dimensions
	plint geoNx = 50, geoNy = 50, geoNz = 50;

	// Number of iterations
	plint Nt = 300000;

	// Number iterations between flow and particle output
	plint output_interval = 100;

	// Input parameters
	double Cas[] = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
	double tau = 1; 			// LBM relaxation parameter
	double Re_p = 0.1;			// Particle Reynolds number (= (rbc radius)^2 * gamma / nu)
	double lambda = 1;			// Viscosity ratio between the fluid inside and outside the rbcs
	//double Ca = 0.1;			// Capillary number (mu*(rbc radius)*gamma / (rbc shear modulus)
	double k_bend_rel = 0.006;	// Relative bending energy (General rule: buckling occurs if > 0.01)
	double poisson_ratio = 0.99;// Poisson ratio (controls the area expansion stiffness, set to ~<1 for area incompressibility)

	// Load mesh (defines the equilibrium shape of the particle)
	ParticleShapeLibrary<double> shape_library;
	shape_library.read_and_store_mesh("sphere_TRImesh_626n_1248el.msh", "RBC");
	ParticleShape<double> * shape = shape_library.get_by_tag("RBC");
	double a_lb = std::sqrt(shape->get_area() / (4. * M_PI));

	// Deduce parameters
	double nu_lb = (tau - 0.5) / 3.0;
	double gamma_lb = Re_p * nu_lb / (a_lb*a_lb);
	double uMax = gamma_lb * geoNy / 2.0;

	// Relaxation frequencies
	double tau_inner = 0.5 + lambda*(tau - 0.5);
	double omega = 1. / tau;
	double omega_inner = 1. / tau_inner;

	for(int i = 0; i < 11; ++i) {
		double Ca = Cas[i];
		std::stringstream ss;
		ss << "rbcInShearCa" << Ca << "_lambda" << lambda << "/";
		std::string folder = ss.str();
		io::mkdir(folder.c_str());
		global::directories().setOutputDir(folder);
		global::directories().setInputDir(global::directories().getOutputDir());

		// Membrane parameters
		double G_lb = a_lb * gamma_lb * nu_lb / Ca;
		double K_area_desired = (poisson_ratio + 1)/(1 - poisson_ratio) * G_lb;
		double K_bend = k_bend_rel * G_lb * shape->get_area();
		RBCParameters<double> rbc_params = create_rbc_params(shape, G_lb, K_area_desired, K_bend);
		RBCParticle<double> rbc(shape_library.get_by_tag("RBC"), rbc_params);

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
		pcout << "  uMax:" << uMax << std::endl;
		pcout << "  tau (outer fluid):" << tau << std::endl;
		pcout << "  tau (inner fluid):" << tau_inner << std::endl;
		pcout << "==============" << std::endl;

		// Create equilbrium shape (only for RBCs)
		shrink_rbc_volume(rbc, (double) 0.59*rbc.shape()->get_volume(), 40000);

		// Rotate it so that its minor axis is parallel to the x-axis
		rbc.set_minor_axis_orientation(Array<double, 3>(1, 0, 0));

		Array<double, 3> channel_mid(0.5*geoNx, 0.5*geoNy, 0.5*geoNz);

		// Create lattice
		MultiBlockLattice3D<double, DESCRIPTOR> lattice(
				geoNx+1, geoNy+1, geoNz+1,
				new GuoExternalForceVOFDynamics<double, DESCRIPTOR>(omega, omega_inner) );
		Box3D domain = lattice.getBoundingBox();
		lattice.toggleInternalStatistics(false);

		OnLatticeBoundaryCondition3D<double,DESCRIPTOR> * boundaryCondition
				= createInterpBoundaryCondition3D<double,DESCRIPTOR>();

		// Setup boundary conditions
		shearFlowSetup(lattice, uMax, *boundaryCondition);

		// Immersed boundary module
		ImmersedBoundaryDynamics3D<double, DESCRIPTOR, Periodicity> fsi(
				lattice,
				shape_library );

		// Velocity and concentration fields
		std::auto_ptr<MultiTensorField3D<double, 3> > velocity = generateMultiTensorField<double, 3>(domain);
		std::auto_ptr<MultiScalarField3D<double> > H = generateMultiScalarField<double>(domain);

		// Add rbc
		// Move RBC to channel center
		//rbc.set_center_of_mass(channel_mid);
		rbc.set_center_of_mass(Array<double, 3>(channel_mid[0], channel_mid[1], channel_mid[2]));
		fsi.add_particle(&rbc);

		// Initialize fsi
		fsi.init();

		pcout << "===============" << std::endl;
		pcout << "Starting coupled ib-lbm iterations" << std::endl;
		for(plint it = 0; it <= Nt; ++it) {
			if((it % 1000) == 0) {
				pcout << "Iteration " << it << " / " << Nt << std::endl;
				for(plint i=0; i < fsi.count_particles(); ++i){
					const ParticleBase3D<double> * p;
					if(fsi.get_particle(i, p))
						dynamic_cast<const RBCParticle<double> *>(p)->print_energies(std::cout);
				}
			}

			computeVelocity(lattice, *velocity, domain);
			setExternalVector(lattice, domain, DESCRIPTOR<double>::ExternalField::forceBeginsAt, Array<double, 3>(0,0,0));
			applyProcessingFunctional(wrap_ibm_dynamics3D(fsi), domain, lattice, *velocity);

			if((it % output_interval) == 0) {
				fsi.write_lightweight_particle_data(it);
			}

			if((it % 1000) == 0) {
				fsi.write_particles_as_vtk(it);

				// Get concentration field
				applyProcessingFunctional(new ComputeConcentrationFunctional<double, DESCRIPTOR>, domain, lattice, *H);

				VtkImageOutput3D<double> vtkOut(createFileName("vtk", it, 6), 1.);
				vtkOut.writeData<3,float>(*velocity, "velocity", 1.);
				vtkOut.writeData<float>(*H, "concentration", 1.);
			}
			lattice.collideAndStream();
		}

		//Profile::write_report("profile");

		parallelIO::save(lattice, "checkpoint_lbm", false);
		fsi.save_checkpoint("checkpoint_fsi");
	}
}

