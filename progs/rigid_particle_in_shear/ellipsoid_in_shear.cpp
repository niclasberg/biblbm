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
#include "../common/calibrate_rbc_shape.h"

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
	plint geoNx = 32, geoNy = 32, geoNz = 32;
	Array<double, 3> channel_mid(0.5*geoNx, 0.5*geoNy, 0.5*geoNz);

	// Number of iterations
	plint Nt = 60001;

	// Number iterations between flow and particle output
	plint output_interval = 50;

	// Input parameters
	double tau = 1; 			// LBM relaxation parameter
	double Re_p = 0.1;			// Particle Reynolds number (= (rbc radius)^2 * gamma / nu)
	double lambda = 1;			// Viscosity ratio between the fluid inside and outside the rbcs
	double Ca = 0.01;			// Capillary number (mu*(rbc radius)*gamma / (rbc shear modulus)
	double k_bend_rel = 0.01;	// Relative bending energy (General rule: buckling occurs if > 0.01)
	double poisson_ratio = 0.9;// Poisson ratio (controls the area expansion stiffness, set to ~1 for area incompressibility)

	// Mesh files
	std::vector<std::pair<std::string, std::string> > mesh_files;
	mesh_files.push_back(std::make_pair("ellipsoid_1.5_3_3_TRImesh_74n_144el.msh", "74N"));
	mesh_files.push_back(std::make_pair("ellipsoid_1.5_3_3_TRImesh_166n_328el.msh", "166N"));
	mesh_files.push_back(std::make_pair("ellipsoid_1.5_3_3_TRImesh_374n_744el.msh", "374N"));

	for(std::vector<std::pair<std::string, std::string> >::iterator mIt = mesh_files.begin(); mIt != mesh_files.end(); ++mIt) {
		// Create output directory
		std::stringstream ss;
		ss << "rigidParticleInShearCa" << Ca << "_Rep" << Re_p << "_N" << mIt->second << "/";
		global::directories().setOutputDir(ss.str());
		global::directories().setInputDir(global::directories().getOutputDir());
		io::mkdir(global::directories().getOutputDir().c_str());

		// Load mesh (defines the equilibrium shape of the particle)
		ParticleShapeLibrary<double> shape_library;
		shape_library.read_and_store_mesh(mIt->first.c_str(), "RBC");
		ParticleShape<double> * shape = shape_library.get_by_tag("RBC");
		double a_lb = std::sqrt(shape->get_area() / (4*M_PI));

		// Deduce parameters
		double nu_lb = (tau - 0.5) / 3.0;
		double gamma_lb = Re_p *nu_lb / (a_lb*a_lb);
		double uMax = gamma_lb * geoNx / 2.0;

		// Relaxation frequencies
		double tau_inner = 0.5 + lambda*(tau - 0.5);
		double omega = 1. / tau;
		double omega_inner = 1. / tau_inner;

		// Membrane parameters
		//SemiRigidParticleParams<double> particle_params;
		
		double G_lb = a_lb * gamma_lb * nu_lb / Ca;
		double K_area_desired = (poisson_ratio + 1)/(1 - poisson_ratio) * G_lb;
		double K_bend = k_bend_rel * G_lb * shape->get_area();

		// In plane properties
		// Get average link length
		/*double lavg = 0;
		for(ParticleShape<double>::link_const_iterator it = shape->links_begin(); it != shape->links_end(); ++it)
			lavg += it->length;
		lavg /= shape->count_links();
		particle_params.l0 = lavg;

		// The shear modulus is proportional to k_in_plane
		particle_params.k_in_plane = 1;
		particle_params.k_in_plane = G_lb / particle_params.shear_modulus();
		particle_params.k_out_of_plane = particle_params.k_in_plane;

		particle_params.k_bend = k_bend_rel * G_lb * shape->get_area();

		SemiRigidParticle3D<double> rbc(shape_library.get_by_tag("RBC"), particle_params);*/

		RBCParticle<double> rbc(shape, create_platelet_params(shape, G_lb, K_area_desired, K_bend));

		pcout << "RBC Parameters:" << std::endl;
		pcout << "  Global area stiffness: " << rbc.params().k_area_global << std::endl;
		pcout << "  Local area stiffness: " << rbc.params().k_area_local << std::endl;
		pcout << "  Global volume stiffness: " << rbc.params().k_volume << std::endl;
		pcout << "  Bending stiffness: " << rbc.params().k_bend << std::endl;
		pcout << "  Equilibrium angle: " << rbc.params().theta0 << std::endl;
		pcout << "  Shear modulus: " << rbc.params().G() << std::endl;
		pcout << "  Area modulus " << rbc.params().K() << std::endl;
		pcout << "  Young's modulus: " << rbc.params().youngs_modulus() << std::endl;
		pcout << "==============" << std::endl;
		pcout << "Non-dimensional parameters:" << std::endl;
		pcout << "  Re_p: " << gamma_lb * util::sqr(a_lb) / nu_lb << std::endl;
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
			if(it.is_multiple_of(500)) {
				for(plint i=0; i < fsi.count_particles(); ++i){
					const ParticleBase3D<double> * p;
					if(fsi.get_particle(i, p)) {
						const RBCParticle<double> * pRBC = dynamic_cast<const RBCParticle<double> *>(p);
						std::cout << "Iteration " << it << " / " << Nt << ", rbc area: " << pRBC->area() << ", volume: " << pRBC->volume() << std::endl;
						pRBC->print_energies(std::cout);
					}
				}
				//pcout << "Iteration " << it << " / " << Nt << " total volume: " << fsi.compute_total_volume() << std::endl;
			}

			computeVelocity(lattice, *velocity, domain);

			if(false && it.is_multiple_of(10)) {
				pcout << "Writing lightweight particle data at iteration " << it << std::endl;
				fsi.write_lightweight_particle_data(it.to_plint());
			}

			if(it.is_multiple_of(500)) {
				pcout << "Writing particle vtks at iteration " << it << std::endl;
				fsi.write_particles_as_vtk(it.to_plint());
			}

			if(it.is_multiple_of(500)) {
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
}

