// Comment away this line to disable profiling:
#define FSI_PROFILE

#include "palabos3D.h"
#include "palabos3D.hh"
#include "fsi/headers.h"
#include "fsi/headers.hh"
#include "fsi/SemiRigidParticle.h"
#include "fsi/SemiRigidParticle.hh"
#include "fsi/Time.h"
#include "../common/shear_flow_setup.h"

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
	global::directories().setOutputDir("test_ell/");
	global::directories().setInputDir(global::directories().getOutputDir());

	// Create output directory
	io::mkdir(global::directories().getOutputDir().c_str());

	// Domain dimensions
	plint geoNx = 32, geoNy = 32, geoNz = 32;
	Array<double, 3> channel_mid(0.5*geoNx, 0.5*geoNy, 0.5*geoNz);

	// Number of iterations
	plint Nt = 30000;

	// Number iterations between flow and particle output
	plint output_interval = 50;

	// Input parameters
	double tau = 1; 			// LBM relaxation parameter
	double Re_p = 0.01;			// Particle Reynolds number (= (rbc radius)^2 * gamma / nu)
	double lambda = 1;			// Viscosity ratio between the fluid inside and outside the rbcs
	double Ca = 0.01;			// Capillary number (mu*(rbc radius)*gamma / (rbc shear modulus)
	double k_bend_rel = 0.001;	// Relative bending energy (General rule: buckling occurs if > 0.01)
	double poisson_ratio = 0.990;// Poisson ratio (controls the area expansion stiffness, set to ~1 for area incompressibility)

	// Load mesh (defines the equilibrium shape of the particle)
	ParticleShapeLibrary<double> shape_library;
	shape_library.read_and_store_mesh("ellipsoid_2_4_4_TRImesh_174n_344el.msh", "RBC");
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
	SemiRigidParticleParams<double> particle_params;
	double G_lb = a_lb * gamma_lb * nu_lb / Ca;

	// In plane properties
	// Get average link length
	double lavg = 0;
	for(ParticleShape<double>::link_const_iterator it = shape->links_begin(); it != shape->links_end(); ++it)
		lavg += it->length;
	lavg /= shape->count_links();
	particle_params.l0 = lavg;

	// The shear modulus is proportional to k_in_plane
	particle_params.k_in_plane = 1;
	particle_params.k_in_plane = G_lb / particle_params.shear_modulus();
	particle_params.k_out_of_plane = particle_params.k_in_plane;

	particle_params.k_bend = k_bend_rel * G_lb * shape->get_area();

	SemiRigidParticle3D<double> rbc(shape_library.get_by_tag("RBC"), particle_params);

	pcout << "RBC Parameters:" << std::endl;
	pcout << "  In-plane shear stiffness: " << particle_params.k_in_plane << std::endl;
	pcout << "  Out of plane shear stiffness: " << particle_params.k_out_of_plane << std::endl;
	pcout << "==============" << std::endl;
	pcout << "Non-dimensional parameters:" << std::endl;
	pcout << "  Re_p: " << gamma_lb * util::sqr(shape->get_radius()) / nu_lb << std::endl;
	pcout << "==============" << std::endl;
	pcout << "Numerics:" << std::endl;
	pcout << "  uMax:" << uMax << std::endl;
	pcout << "  tau (outer fluid):" << tau << std::endl;
	pcout << "  tau (inner fluid):" << tau_inner << std::endl;
	pcout << "==============" << std::endl;

	// Move RBC to channel center
	rbc.set_center_of_mass(channel_mid);

	// Rotate it so that its minor axis is parallel to the x-axis
	//rbc.set_minor_axis_orientation(Array<double, 3>(1, 0, 0));
	//rbc.set_major_axis_orientation(Array<double, 3>(1, 0, 0));

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

	pcout << "===============" << std::endl;
	pcout << "Starting coupled ib-lbm iterations" << std::endl;
	fsi::Time it(0);
	bool first_iteration = true;

	// Main computation loop
	for(; it <= Nt; ++it) {
		if(it.is_multiple_of(100)) {
			pcout << "Iteration " << it << " / " << Nt << " total volume: " << fsi.compute_total_volume() << std::endl;
		}

		computeVelocity(lattice, *velocity, domain);

		/*if(it.is_multiple_of(particle_lw_interval)) {
			pcout << "Writing lightweight particle data at iteration " << it << std::endl;
			fsi.write_lightweight_particle_data(it.to_plint());
		}*/

		if(it.is_multiple_of(1000)) {
			pcout << "Writing particle vtks at iteration " << it << std::endl;
			fsi.write_particles_as_vtk(it.to_plint());
		}

		if(it.is_multiple_of(1000)) {
			pcout << "Writing flow field at iteration " << it << std::endl;

			// Get concentration and force field
			//applyProcessingFunctional(new ComputeConcentrationFunctional<double, DESCRIPTOR>, domain, lattice, *H);

			VtkImageOutput3D<double> vtkOut(createFileName("vtk", it.to_plint(), 6), 1.);
			vtkOut.writeData<3,float>(*velocity, "velocity", 1.);
			//vtkOut.writeData<float>(*H, "concentration", 1.);
		}

		// Perform IB step
		setExternalVector(lattice, domain, DESCRIPTOR<double>::ExternalField::forceBeginsAt, Array<double, 3>(0., 0., 0.));
		applyProcessingFunctional(wrap_ibm_dynamics3D(fsi), domain, lattice, *velocity);

		// Do lattice boltzmann step
		Profile::start_timer("lbm");
		lattice.collideAndStream();
		Profile::stop_timer("lbm");

		// Just a flag
		first_iteration = false;
	}

	Profile::write_report("profile");

	parallelIO::save(lattice, "checkpoint_lbm", false);
	fsi.save_checkpoint("checkpoint_fsi");
}

