// Comment away this line to disable profiling:
#define FSI_PROFILE

#include "palabos3D.h"
#include "palabos3D.hh"
#include "fsi/headers.h"
#include "fsi/headers.hh"
#include "../common/shear_flow_setup.h"

#include <algorithm>

using namespace plb;
using namespace fsi;

#define DESCRIPTOR descriptors::ForcedPhaseD3Q19Descriptor

// Define the periodicity of the domain (periodicity_x, periodicity_y, periodicity_z)
typedef Periodicity3D<double, true, false, true> Periodicity;

int main(int argc, char ** argv)
{
	// Initialize palabos
	plbInit(&argc, &argv);

	// Input parameters default values
	plint N = 120;				// Domain dimensions
	double tau = 0.8;			// LBM relaxation parameter
	double Re_p = 0.1;			// Particle Reynolds number (= (rbc radius)^2 * gamma / nu)
	double lambda = 1;			// Viscosity ratio between the fluid inside and outside the rbcs
	double Ca = 0.15;			// Capillary number (mu*(rbc radius)*gamma / (rbc shear modulus)
	double k_bend_rel = 0;		// Relative bending energy
	double poisson_ratio = 1;	// Poisson ratio (controls the area expansion stiffness, set to ~1 for area incompressibility)

	// Read command line parameters
	{
		std::string param;
		for(int i = 1; i < plb::global::argc(); i+=2) {
			try {
				global::argv(i).read(param);
				if(param.compare("-tau") == 0)
					global::argv(i+1).read(tau);
				else if(param.compare("-Ca") == 0)
					global::argv(i+1).read(Ca);
				else if(param.compare("-Re") == 0)
					global::argv(i+1).read(Re_p);
				else if(param.compare("-N") == 0)
					global::argv(i+1).read(N);
				else if(param.compare("-lambda") == 0)
					global::argv(i+1).read(lambda);
				else
					throw PlbIOException("Unknown parameter");
			} catch(PlbIOException& exception) {
				// Print the corresponding error message.
				pcout << exception.what() << std::endl;
				pcout << "Usage: capsule_in_shear [-tau value] [-Ca value] [-Re value] [-N value] [-lambda value]" << std::endl;
				return EXIT_FAILURE;
			}
		}
	}

	plint geoNx = N, geoNy = N, geoNz = N;
	Array<double, 3> channel_mid(0.5*geoNx, 0.5*geoNy, 0.5*geoNz);

	// Setup input/output
	global::IOpolicy().activateParallelIO(true);
	{
		std::stringstream ss;
		ss << "spheroid_N" << N << "_Ca" << Ca << "_tau" << tau << "/";
		std::string folder = ss.str();
		std::replace(folder.begin(), folder.end(), '.', '_');
		global::directories().setOutputDir(folder);
		global::directories().setInputDir(global::directories().getOutputDir());
	}

	// Create output directory
	io::mkdir(global::directories().getOutputDir().c_str());

	// Load mesh (defines the equilibrium shape of the particle)
	ParticleShapeLibrary<double> shape_library;
	if(N == 60)
		shape_library.read_and_store_mesh("sphere_r6_TRImesh_1046n_2088el.msh", "RBC");
	else if(N == 120)
		shape_library.read_and_store_mesh("sphere_r12_TRImesh_4180n_8356el.msh", "RBC");
	else
		throw std::runtime_error("Unsupported resolution");
	ParticleShape<double> * shape = shape_library.get_by_tag("RBC");
	double a_lb = shape->get_radius();

	// Deduce parameters
	double nu_lb = (tau - 0.5) / 3.0;
	double gamma_lb = Re_p *nu_lb / (a_lb*a_lb);
	double uMax = gamma_lb * geoNx / 2.0;

	// Number of iterations (100 / gamma)
	plint Nt = std::ceil(20 / gamma_lb);

	// Number iterations between flow and particle output
	plint output_interval = std::ceil(1 /gamma_lb / 20);

	// Relaxation frequencies
	double tau_inner = 0.5 + lambda*(tau - 0.5);
	double omega = 1. / tau;
	double omega_inner = 1. / tau_inner;

	// Membrane parameters
	CapsuleParameters<double> rbc_params;
	double G_lb = a_lb * gamma_lb * nu_lb / Ca;

	// Shear modulus
	rbc_params.G = G_lb;

	// Poisson ratio
	rbc_params.nu_s = poisson_ratio;

	// Area compression ratio / shear modulus (only for Skalak model)
	rbc_params.C = 1.;

	// Bending stiffness
	rbc_params.kb = k_bend_rel * G_lb * shape->get_area();

	DeformableCapsuleParticle<double> rbc(shape_library.get_by_tag("RBC"), rbc_params);

	pcout << "RBC Parameters:" << std::endl;
	pcout << "==============" << std::endl;
	pcout << "Non-dimensional parameters:" << std::endl;
	pcout << "  Re_p: " << gamma_lb * util::sqr(shape->get_radius()) / nu_lb << std::endl;
	pcout << "  Ca_G: " << gamma_lb * shape->get_radius() * nu_lb / rbc_params.G << std::endl;
	pcout << "  lambda: " << lambda << std::endl;
	pcout << "  Relative bending stiffness (kb/GA): " << rbc_params.kb / (rbc_params.G * rbc.area()) << std::endl;
	pcout << "  Poisson ratio " << rbc_params.nu_s << std::endl;
	pcout << "==============" << std::endl;
	pcout << "Numerics:" << std::endl;
	pcout << "  uMax:" << uMax << std::endl;
	pcout << "  tau (outer fluid):" << tau << std::endl;
	pcout << "  tau (inner fluid):" << tau_inner << std::endl;
	pcout << "==============" << std::endl;
	pcout << "Domain" << std::endl;
	pcout << "  " << geoNx << " x " << geoNy << " x " << geoNz << std::endl;
	pcout << "  Total number of fluid nodes: " << (geoNx+1)*(geoNy+1)*(geoNz+1) << std::endl;

	// Move RBC to channel center
	rbc.set_center_of_mass(channel_mid);

	// Rotate it so that its minor axis is parallel to the x-axis
	//rbc.set_minor_axis_orientation(Array<double, 3>(1, 0, 0));
	rbc.set_major_axis_orientation(Array<double, 3>(0, 1, 0));

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
	fsi.add_particle(&rbc);

	// Initialize fsi
	fsi.init();

	// Create output file
	std::fstream fout;
	{
		std::stringstream ss;
		ss << global::directories().getOutputDir() << "deformation_p" << global::mpi().getRank() << ".txt";
		fout.open(ss.str().c_str(), std::ios::out);
	}

	pcout << "===============" << std::endl;
	pcout << "Starting coupled ib-lbm iterations" << std::endl;

	for(plint it = 0; it <= Nt; ++it) {
		if((it % 50) == 0)
			pcout << "Iteration " << it << " / " << Nt << ", particle volume = " << fsi.compute_total_volume() << std::endl;

		computeVelocity(lattice, *velocity, domain);

		if((it % output_interval) == 0) {
			fsi.write_particles_as_vtk(it);

			// Get concentration field
			//applyProcessingFunctional(new ComputeConcentrationFunctional<double, DESCRIPTOR>, domain, lattice, *H);

			VtkImageOutput3D<double> vtkOut(createFileName("vtk", it, 6), 1.);
			vtkOut.writeData<3,float>(*velocity, "velocity", 1.);
			//vtkOut.writeData<float>(*H, "concentration", 1.);
		}

		setExternalVector(lattice, domain, DESCRIPTOR<double>::ExternalField::forceBeginsAt, Array<double, 3>(0,0,0));
		applyProcessingFunctional(wrap_ibm_dynamics3D(fsi), domain, lattice, *velocity);

		// Write deformation state
		if((it % 10) == 0) {
			ParticleBase3D<double> * p;
			if(fsi.get_particle(0, p)) {
				// Deformation tensor
				Matrix<double, 3> moment_of_inertia;
				dynamic_cast<DeformableCapsuleParticle<double> *>(p)->compute_moment_of_intertia(moment_of_inertia);

				fout << gamma_lb*it << " " << p->volume();
				for(int i = 0; i < 3; ++i)
					for(int j = 0; j < 3; ++j)
						fout << " " << moment_of_inertia(i, j);
				fout << std::endl;
			}
		}

		// Compute forces and update node positions
		Profile::start_timer("lbm");
		lattice.collideAndStream();
		Profile::stop_timer("lbm");
	}
	fout.close();

	//parallelIO::save(lattice, "checkpoint_lbm", false);
	//fsi.save_checkpoint("checkpoint_fsi");
}

