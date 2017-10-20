//#define FSI_PROFILE

#include "palabos3D.h"
#include "palabos3D.hh"
#include "fsi/headers.h"
#include "fsi/headers.hh"
#include "../common/calibrate_rbc_shape.h"
#include "../common/util.h"
#include "../common/rbc_parameters.h"
#include "fsi/Time.h"
#include "fsi/GuoBoundaryCondition.h"
#include "fsi/GuoBoundaryCondition.hh"
#include "fsi/ParticlePositionInitializer.h"
#include "fsi/ParticlePositionInitializer.hh"

using namespace plb;
using namespace fsi;

#define DESCRIPTOR descriptors::ForcedPhaseD3Q19Descriptor
typedef Periodicity3D<double, true, false, false> Periodicity;

// Margin around the boundary (needed for boundary conditions)
plint margin = 1;
// The parallel envelope of the fluid (2 is needed for the off-lattice boundary condition).
plint extendedEnvelopeWidth = 2;

// Domain
plint diameter;						// Diameter in lattice units
plint length; 						// Pipe length

bool load_checkpoint;
plint Nt;							// Maximum number of iterations
plint flow_vtk_interval;			// Iterations between output
plint particle_vtk_interval;
plint particle_lw_interval;
plint checkpoint_interval;			// Iterations between checkpoint output
double omega, omega_inner; 			// Relaxation frequencies'
double uMax;
RBCParameters<double> rbc_params;	// Red blood cell membrane parameters
SemiRigidParticleParams<double> platelet_params;
Array<double, 3> external_force; 	// Pressure gradient
PipeBoundary<double> pipe_boundary;
bool second_order_bcs = true;		// true  => use Guo extrapolation boundary conditions for cylinder surface
									// false => use staircase bounce back at cylinder surface
double rbc_conc;
double platelet_conc;
std::string output_folder;
std::string checkpoint_folder;
PPInteractionForce<double> * pp_interaction = 0;
ForceDecorator<double> * wall_interaction = 0;

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
		T r_sqr = (plb::util::sqr(iY - this->boundary.y0) + plb::util::sqr(iZ - this->boundary.z0)) / plb::util::sqr(this->boundary.radius);
		//u[0] = r_sqr < plb::util::sqr(boundary.radius) ? (T) lattice_u * ((T) 1 - r_sqr) : (T)0.;
		u[0] = (T) 0.;
		u[1] = (T) 0;
		u[2] = (T) 0;
	}
private:
	const plb::fsi::PipeBoundary<T> & boundary;
	T lattice_u;
};

template<class T>
class BoundarySurface : public plb::DomainFunctional3D {
public:
	BoundarySurface(const Boundary<T> & boundary_) : boundary(boundary_) { }

	virtual BoundarySurface<T> * clone() const { return new BoundarySurface(*this); }

	virtual bool operator()(plint iX, plint iY, plint iZ) const
	{
		return ! boundary.contains(Array<T, 3>(iX, iY, iZ));
	}

private:
	const Boundary<T> & boundary;
};

template<class T>
MorsePotential<T> create_morse_potential(const XMLreaderProxy & node)
{
	double r_min, beta, strength;
	node["r_min"].read(r_min);
	node["beta"].read(beta);
	node["strength"].read(strength);
	return MorsePotential<T>(r_min, beta, strength);
}

template<class T>
SpringPotential<T> create_spring_potential(const XMLreaderProxy & node)
{
	double r_min, strength;
	node["r_min"].read(r_min);
	node["strength"].read(strength);
	return SpringPotential<T>(r_min, strength);
}

void read_parameters(ParticleShapeLibrary<double> & shape_library, std::string file_name)
{
	// Input parameters
	double tau; 			// Relaxation parameter
	double Re_d;			// Channel Reynolds number
	double lambda;			// Viscosity ratio between the fluid inside and outside the rbcs
	double Ca;				// Capillary number (mu*(rbc radius)*gamma / (rbc shear modulus)
	double Ca_p;			// Capillary number (mu*(platelet radius)*gamma / (platelet shear modulus)
	double k_bend_rel;		// Relative bending energy
	double poisson_ratio;	// Poisson ratio
	std::string rbc_mesh;
	std::string platelet_mesh;
	std::string wall_interaction_type, pp_interaction_type;

	try {
		// Open the XMl file
		XMLreader xmlFile(file_name.c_str());

		// Pipe geometry
		xmlFile["geometry"]["diameter"].read(diameter);
		xmlFile["geometry"]["length"].read(length);

		// Set boundary dimensions
		pipe_boundary.y0 = 0.;
		pipe_boundary.z0 = 0.;
		pipe_boundary.radius = (double) diameter / 2.;
		pipe_boundary.length = (double) length;

		// Reynolds number
		xmlFile["physics"]["Re_d"].read(Re_d);

		// Read rbc properties
		xmlFile["particles"]["rbc"]["Ca"].read(Ca);
		xmlFile["particles"]["rbc"]["mesh"].read(rbc_mesh);
		xmlFile["particles"]["rbc"]["poisson_ratio"].read(poisson_ratio);
		xmlFile["particles"]["rbc"]["Xi_b"].read(k_bend_rel);
		xmlFile["particles"]["rbc"]["lambda"].read(lambda);
		xmlFile["particles"]["rbc"]["concentration"].read(rbc_conc);

		// Read platelet properties
		xmlFile["particles"]["platelet"]["concentration"].read(platelet_conc);
		xmlFile["particles"]["platelet"]["mesh"].read(platelet_mesh);
		xmlFile["particles"]["platelet"]["Ca"].read(Ca_p);

		// Read wall interaction parameters
		xmlFile["particles"]["wall_interaction"]["type"].read(wall_interaction_type);
		if(wall_interaction_type.compare("morse") == 0) {
			wall_interaction =  new WallInteraction<double, MorsePotential<double> >(
					pipe_boundary,
					create_morse_potential<double>(xmlFile["particles"]["wall_interaction"]));
		} else if(wall_interaction_type.compare("spring") == 0) {
			wall_interaction =  new WallInteraction<double, SpringPotential<double> >(
					pipe_boundary,
					create_spring_potential<double>(xmlFile["particles"]["wall_interaction"]));
		} else if(wall_interaction_type.compare("none") != 0) {
			std::cerr << "Unknown wall interaction type " << wall_interaction_type << std::endl;
		}

		// Read particle-particle interaction parameters
		xmlFile["particles"]["particle_interaction"]["type"].read(pp_interaction_type);
		if(pp_interaction_type.compare("morse") == 0) {
			pp_interaction =  new PPPotentialForce<double, MorsePotential<double> >(
					create_morse_potential<double>(xmlFile["particles"]["particle_interaction"]));
		} else if(pp_interaction_type.compare("spring") == 0) {
			pp_interaction =  new PPPotentialForce<double, SpringPotential<double> >(
					create_spring_potential<double>(xmlFile["particles"]["particle_interaction"]));
		} else if(pp_interaction_type.compare("none") != 0) {
			std::cerr << "Unknown particle interaction type " << pp_interaction_type << std::endl;
		}

		// Read numerics
		xmlFile["numerics"]["tau"].read(tau);
		xmlFile["numerics"]["max_iterations"].read(Nt);

		// Checkpoints
		xmlFile["checkpoints"]["folder"].read(checkpoint_folder);
		xmlFile["checkpoints"]["start_from_checkpoint"].read(load_checkpoint);
		xmlFile["checkpoints"]["write_interval"].read(checkpoint_interval);

		// Append a trailing slash
		if(checkpoint_folder[checkpoint_folder.size()-1] != '/')
			checkpoint_folder.push_back('/');

		// Output
		xmlFile["output"]["folder"].read(output_folder);
		xmlFile["output"]["flow_vtk_interval"].read(flow_vtk_interval);
		xmlFile["output"]["particle_vtk_interval"].read(particle_vtk_interval);
		xmlFile["output"]["particle_lw_interval"].read(particle_lw_interval);

		// Append a trailing slash
		if(output_folder[output_folder.size()-1] != '/')
			output_folder.push_back('/');

	} catch (PlbIOException& exception) {
		pcout << exception.what() << std::endl;
		exit(-1);
	}

	// Compute fluid viscosity and relaxation frequencies
	double nu_lb = (tau - 0.5) / 3.0;
	double tau_inner = 0.5 + lambda*(tau - 0.5);
	omega = 1. / tau;
	omega_inner = 1. / tau_inner;

	// Compute maximal velocity ( = 2 * u_mean)
	uMax = 2. * nu_lb * Re_d / (double) diameter;

	// Load platelet mesh
	shape_library.read_and_store_mesh(platelet_mesh.c_str(), "PLATELET");
	ParticleShape<double> * platelet_shape = shape_library.get_by_tag("PLATELET");	

	// Load RBC mesh
	shape_library.read_and_store_mesh(rbc_mesh.c_str(), "RBC");
	ParticleShape<double> * shape = shape_library.get_by_tag("RBC");
	double a_lb = std::sqrt(shape->get_area() / (4. * M_PI));

	// Deduce parameters
	double gamma_lb = 4. * uMax / (double) diameter;

	// Membrane parameters
	double G_lb = a_lb * gamma_lb * nu_lb / Ca;
	double K_area = (poisson_ratio + 1)/(1 - poisson_ratio) * G_lb;
	double K_bend = k_bend_rel * G_lb * shape->get_area();
	rbc_params = create_rbc_params(shape, G_lb, K_area, K_bend);

	// Create platelet params
	{
		// Compute average link length
		double lavg = 0;		
		for(ParticleShape<double>::link_const_iterator it = platelet_shape->links_begin(); 
				it != platelet_shape->links_end(); ++it)
			lavg += it->length;
		lavg /= platelet_shape->count_links();
		platelet_params.l0 = lavg;
		
		// Compute shear modulus
		double a_p_lb = std::sqrt(platelet_shape->get_area() / (4. * M_PI));
		double G_plat_lb = a_p_lb * gamma_lb * nu_lb / Ca_p;

		// Set in plane energy (proportional to the shear modulus)
		platelet_params.k_in_plane = 1.;
		platelet_params.k_in_plane = G_lb / platelet_params.shear_modulus();
		platelet_params.k_out_of_plane = 10*platelet_params.k_in_plane;
	}

	// External force
	external_force[0] = 4. * nu_lb * gamma_lb / (double) diameter;
	external_force[1] = 0.;
	external_force[2] = 0.;

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
	pcout << "  Re_d: " << uMax * (double) diameter / (2.*nu_lb) << std::endl;
	pcout << "  Re_p: " << gamma_lb * util::sqr(shape->get_radius()) / nu_lb << std::endl;
	pcout << "  Ca_G: " << gamma_lb * shape->get_radius() * nu_lb / rbc_params.G() << std::endl;
	pcout << "  lambda: " << lambda << std::endl;
	pcout << "  Relative bending stiffness (kb/GA): " << rbc_params.k_bend / (rbc_params.G() * shape->get_area()) << std::endl;
	pcout << "  Poisson ratio " << rbc_params.poisson_ratio() << std::endl;
	pcout << "==============" << std::endl;
	pcout << "Numerics:" << std::endl;
	pcout << "  uMax:" << uMax << std::endl;
	pcout << "  tau (outer fluid):" << tau << std::endl;
	pcout << "  tau (inner fluid):" << tau_inner << std::endl;
	pcout << "==============" << std::endl;
}

template<class T>
void setup_flow(
		MultiBlockLattice3D<T,DESCRIPTOR> *& fluidLattice,
		GuoRigidWallBoundary<T, DESCRIPTOR> *& boundary_condition)
{
	 // Allocate the data for the fluid.
	fluidLattice = generateMultiBlockLattice<T,DESCRIPTOR> (
			pipe_boundary.get_bounding_box(margin),
			new GuoExternalForceVOFDynamics<T, DESCRIPTOR>(omega, omega_inner),
			extendedEnvelopeWidth).release();

	// Apply periodic boundary condition in the x-direction
	fluidLattice->periodicity().toggle(0, Periodicity::get_x());
	fluidLattice->periodicity().toggle(1, Periodicity::get_y());
	fluidLattice->periodicity().toggle(2, Periodicity::get_z());
	fluidLattice->toggleInternalStatistics(false);

	// Insert cylinder wall boundary condition
	if(second_order_bcs) {
		boundary_condition = new GuoRigidWallBoundary<T, DESCRIPTOR>(pipe_boundary, *fluidLattice);
		boundary_condition->insert();
	} else {
		defineDynamics<T,DESCRIPTOR>(*fluidLattice,
				fluidLattice->getBoundingBox(),
				new BoundarySurface<T>(pipe_boundary),
				new BounceBack<T,DESCRIPTOR>() );
	}

	// Optimization: turn off dynamics for the cells outside of the cylindrical domain
	defineDynamics<T,DESCRIPTOR>(
				*fluidLattice,
				fluidLattice->getBoundingBox(),
				new BoundaryOutside<T>(pipe_boundary, std::sqrt(2)),
				new plb::NoDynamics<T,DESCRIPTOR>() );

	// Read checkpoint file if requested
	if(load_checkpoint)
		parallelIO::load("checkpoint_lbm", *fluidLattice, false);
	else {
		setExternalVector(*fluidLattice,
				fluidLattice->getBoundingBox(),
				DESCRIPTOR<T>::ExternalField::forceBeginsAt,
				external_force);

		// Default-initialize at equilibrium with zero-velocity and fixed pressure.
		initializeAtEquilibrium(*fluidLattice,
				fluidLattice->getBoundingBox(),
				PipeFlowDensityAndVelocity<double>(pipe_boundary, uMax));

		// Export boundary geometry
		pipe_boundary.writeVTK(global::directories().getOutputDir() + std::string("boundary.vtu"));
	}

	// Execute all data processors once to start the simulation off with well-defined initial values.
	fluidLattice->initialize();
}

template<class T>
void create_particles_regular(
		ImmersedBoundaryDynamics3D<T, DESCRIPTOR, Periodicity> & fsi,
		MultiBlockLattice3D<T,DESCRIPTOR> & lattice,
		const ParticleShapeLibrary<T> & shape_library)
{
	// Create intial condition for particles
	if(load_checkpoint) {
		pcout << "Loading particle checkpoint" << std::endl;
		fsi.load_checkpoint("checkpoint_fsi");
	} else {
		RBCParticle<double> rbc(shape_library.get_by_tag("RBC"), rbc_params);

		// Create equilbrium shape
		shrink_rbc_volume(rbc, (double) 0.59*rbc.shape()->get_volume(), 10000);

		rbc.set_minor_axis_orientation(Array<T, 3>(1, 0, 0));
		rbc.set_center_of_mass(Array<T, 3>(0, 0, 0));

		pcout << "Creating initial particle distribution" << std::endl;

		// Compute the number of particles (concentration*domain volume / rbc_volume)
		plint num_particles = std::ceil((rbc_conc * M_PI * util::sqr(pipe_boundary.radius) * pipe_boundary.length) / rbc.volume());
		std::vector<Array<double, 3> > positions;

		if(global::mpi().isMainProcessor()) {
			double size_x = rbc.bounding_box().x1 - rbc.bounding_box().x0;
			double size_r = std::max(rbc.bounding_box().y1 - rbc.bounding_box().y0,
									 rbc.bounding_box().z1 - rbc.bounding_box().z0);

			// Number of particles that fits within one radius (minus a margin of 0.5 so that fsi will still work)
			double fsi_margin = 0.5;
			plint nr = std::floor((pipe_boundary.radius + 0.5*size_r - fsi_margin) / size_r);
			double unfilled_r = pipe_boundary.radius - fsi_margin - (nr-0.5)*size_r;
			double spacing_r = unfilled_r / (double) nr;		//empty radial space between each particle

			if(nr <= 0) {
				std::cerr << "The domain is too small, particles cannot fit. Pipe radius = " << pipe_boundary.radius << ", particle diameter = " << size_r << std::endl;
				exit(-1);
			}

			// Create radial positions
			std::vector<double> rs(nr);
			rs[0] = 0;
			for(plint i = 1; i < nr; ++i)
				rs[i] = rs[i-1] + size_r + spacing_r;

			// Determine how many particles that can be fitted at each radial position
			std::vector<plint> nthetas(nr);
			std::vector<double> theta_spacing(nr);
			std::vector<double> dthetas(nr);
			nthetas[0] = 1;
			theta_spacing[0] = 0;
			dthetas[0] = 0;
			for(plint i = 1; i < nr; ++i) {
				// Some simple trigonometry to determine how much angular space is occupied by each rbc
				double r = rs[i] / (size_r/2.);
				double particle_dtheta = 2. * std::atan2(std::sqrt(1-(1./r/r)), r - 1./r);

				nthetas[i] = std::floor(2*M_PI / particle_dtheta);
				theta_spacing[i] = (2.*M_PI - particle_dtheta*nthetas[i]) / (double) (nthetas[i]);
				dthetas[i] = particle_dtheta + theta_spacing[i];
			}

			// Compute total number of particles / x-position
			plint particles_per_x = std::accumulate(nthetas.begin(), nthetas.end(), 0);

			// Determine how many x-positions are required
			plint nx = num_particles / particles_per_x; //integer division, rounds towards 0
			if(particles_per_x*nx < num_particles)
				++nx;

			double x_spacing = (pipe_boundary.length - nx*size_x) / (double) nx;

			if(x_spacing < 0) {
				std::cerr << "Could not generate an intial distribution, the rbc concentration is too high. " << std::endl;
				exit(-1);
			}

			double delta_x = size_x + x_spacing;

			// Distribute the particles
			plint particles_put = 0;
			for(plint ix = 0; ix < nx; ++ix) {
				double x0 = ix*delta_x;
				for(plint ir = 0; ir < nr; ++ir) {
					for(plint itheta = 0; itheta < nthetas[ir]; ++itheta) {
						// Add random perturbations to theta and x
						double theta0 = itheta * dthetas[ir];

						double rands[3];
						rands[0] = rand() / (double) RAND_MAX;
						rands[1] = rand() / (double) RAND_MAX;
						rands[2] = rand() / (double) RAND_MAX;

						double theta = theta0;;
						double x = x0 + 0.9*x_spacing*(rands[1] - 0.5);
						double r = rs[ir] + 0.45*rands[2] * spacing_r;

						// Add to domain
						positions.push_back(Array<double, 3>(x, rs[ir]*std::cos(theta), rs[ir]*std::sin(theta)));

						++particles_put;
					}
				}
			}

			// Randomly remove exessive particles
			while(positions.size() > num_particles) {
				plint ind = rand() % (positions.size() - 1);
				positions.erase(positions.begin()+ind);
			}
		} else {
			positions.resize(num_particles);
		}

		global::mpi().bCast(reinterpret_cast<T*>(&(positions[0])), 3*num_particles);

		for(plint i = 0; i < num_particles; ++i) {
			rbc.set_center_of_mass(positions[i]);
			fsi.add_particle(&rbc);
		}
	}
}

template<class T>
void create_particles(
		ImmersedBoundaryDynamics3D<T, DESCRIPTOR, Periodicity> & fsi,
		MultiBlockLattice3D<T,DESCRIPTOR> & lattice,
		ParticleShapeLibrary<T> & shape_library)
{
	// Create intial condition for particles
	if(load_checkpoint) {
		pcout << "Loading particle checkpoint" << std::endl;
		fsi.load_checkpoint("checkpoint_fsi");
	} else {
		plint Nit = 40000;
		T initial_scale = 0.2;

		// Create RBC and platelet
		RBCParticle<T> rbc(shape_library.get_by_tag("RBC"), rbc_params);
		SemiRigidParticle3D<T> platelet(shape_library.get_by_tag("PLATELET"), platelet_params);

		// Create equilbrium shape
		shrink_rbc_volume(rbc, (double) 0.59*rbc.shape()->get_volume(), 100000);
		rbc.set_minor_axis_orientation(Array<T, 3>(0, 0, 1));
		rbc.set_center_of_mass(Array<T, 3>(0, 0, 0));

		pcout << "Generating initial distribution" << std::endl;

		// Compute the number of particles (concentration*domain volume / rbc_volume)
		T volume = M_PI * util::sqr(pipe_boundary.radius) * pipe_boundary.length;
		plint num_rbcs = std::ceil(rbc_conc *  volume / rbc.volume());
		plint num_platelets = std::ceil(platelet_conc * volume / platelet.volume());
		plint num_particles = num_rbcs + num_platelets;

		// Store the equilibrium shape
		rbc.store_shape(shape_library, "RBC_RIGID");
		ParticleShape<T> * rbc_rigid_shape = shape_library.get_by_tag("RBC_RIGID");

		// Create a new immersed boundary dynamics object for the rigid body simulation
		ImmersedBoundaryDynamics3D<T, DESCRIPTOR, Periodicity> fsi2(lattice, shape_library);

		// Create a boundary that is a little bit smaller than the actual one,
		// this ensures that all particles are within the domain
		PipeBoundary<T> boundary2 = pipe_boundary;
		boundary2.radius = pipe_boundary.radius - 1.5;

		// Create a rigid particle RBC
		RigidParticle3D<T> rbc_rigid(rbc_rigid_shape);
		rbc_rigid.scale() = initial_scale;
		rbc_rigid.density() = 1.;
		rbc_rigid.damping() = 0.1;
		rbc_rigid.update();
	
		// Create a rigid platelet
		RigidParticle3D<T> platelet_rigid(shape_library.get_by_tag("PLATELET"));
		platelet_rigid.scale() = initial_scale;
		platelet_rigid.density() = 1.;
		platelet_rigid.damping() = 0.1;
		platelet_rigid.update();

		// Generate initial distribution
		std::vector<Array<T, 3> > positions;
		std::vector<Quaternion<T> > orientations;
		if(global::mpi().isMainProcessor()) {
			// Generate positions
			double rbc_radius = initial_scale*rbc_rigid_shape->get_radius() + 1.;
			double platelet_radius = initial_scale*shape_library.get_by_tag("PLATELET")->get_radius() + 1.;
			ParticlePositionInitializer<T> position_initializer(
					lattice.getBoundingBox(),
					&boundary2, rbc_radius);
			position_initializer.generate_points(num_platelets, platelet_radius, 'P');
			position_initializer.generate_points(num_rbcs, rbc_radius, 'R');
			
			// Generate positions
			std::vector<Array<T, 3> > platelet_positions;
			position_initializer.get_points('R', positions);
			position_initializer.get_points('P', platelet_positions);
			positions.insert(positions.end(), platelet_positions.begin(), platelet_positions.end());

			// Generate orientations
			for(plint i = 0; i < num_particles; ++i) {
				Array<T, 3> dir((T)std::rand() / RAND_MAX,
						(T)std::rand() / RAND_MAX,
						(T)std::rand() / RAND_MAX);
				T dir_norm = norm(dir);
				if(dir_norm > std::numeric_limits<T>::epsilon()) {
					T theta = 2. * M_PI * (T)std::rand() / RAND_MAX;
					orientations.push_back(fsi::Quaternion<T>(theta, dir/dir_norm));
				} else {
					orientations.push_back(fsi::Quaternion<T>());
				}
			}
		} else {
			positions.resize(num_particles);
			orientations.resize(num_particles);
		}

		// Broadcast to all processors
		global::mpi().bCast(&initial_scale, 1);
		global::mpi().bCast(reinterpret_cast<T*>(&(positions[0])), num_particles*3);
		global::mpi().bCast(reinterpret_cast<T*>(&(orientations[0])), num_particles*4);

		pcout << "Initial distribution created with the RBCs at " << 100*initial_scale << "% of their full size" << std::endl;
		pcout << "Growing RBCs to their full size" << std::endl;

		// Create and insert particles to the fsi object
		for(plint i = 0; i < num_rbcs; ++i) {
			rbc_rigid.scale() = initial_scale;
			rbc_rigid.set_center_of_mass(positions[i]);
			rbc_rigid.orientation() = orientations[i];
			rbc_rigid.update();
			//pcout << orientations[i] << std::endl;
			//pcout << rbc_rigid.get_node(0).pos[0] << std::endl;
			fsi2.add_particle(&rbc_rigid);
		}

		for(plint i = num_rbcs; i < num_particles; ++i) {
			platelet_rigid.set_center_of_mass(positions[i]);
			platelet_rigid.orientation() = orientations[i];
			platelet_rigid.update();
			fsi2.add_particle(&platelet_rigid);
		}

		// Initialize fsi
		fsi2.init();

		// Potential force
		SpringPotential<T> potential(1., 5e-2);

		// Grow the particles while under the influence of collision forces
		T scale;
		pipe_boundary.writeVTK(global::directories().getOutputDir() + std::string("boundary.vtu"));
		for(plint it = 0; it <= Nit; ++it) {
			if((it % 100) == 0)
				pcout << "  Iteration " << it << " / " << Nit << std::endl;
			//if((it % 500) == 0) fsi2.write_particles_as_vtk(it);

			if((it % 10) == 0) {
				// Rescale all particles
				T fraction = (T)it / (T) Nit;
				scale = initial_scale + (1. - initial_scale)*fraction;

				for(typename ImmersedBoundaryDynamics3D<T, DESCRIPTOR, Periodicity>::ObjMapIterator iter = fsi2.particles_begin();
						iter != fsi2.particles_end(); ++iter) {
					plb::fsi::RigidParticle3D<T> * p = dynamic_cast<plb::fsi::RigidParticle3D<T> *>(iter->second);
					if(p) {
						p->scale() = scale;
						p->update();
					}
				}

				fsi2.synchronize_particle_states();
			}

			fsi2.compute_collision_forces(potential, &boundary2);
			fsi2.move_vertices();
		}

		// Scatter all particle positions and orientations
		CommunicationBuffer comm_buff;
		std::vector<plint> proc_list;
		for(plint i = 0; i < global::mpi().getSize(); ++i)
			if(i != global::mpi().getRank())
				proc_list.push_back(i);
		comm_buff.set_proc_list(proc_list);

		for(ImmersedBoundaryDynamics3D<double, DESCRIPTOR, Periodicity>::ObjMapIterator iter = fsi2.particles_begin();
						iter != fsi2.particles_end(); ++iter) {
			plb::fsi::RigidParticle3D<T> * p = dynamic_cast<plb::fsi::RigidParticle3D<T> *>(iter->second);
			for(plint i = 0; i < proc_list.size(); ++i) {
				comm_buff.pack<plint>(proc_list[i], p->get_id());
				comm_buff.pack(proc_list[i], p->center_of_mass());
				comm_buff.pack(proc_list[i], p->orientation());
			}

			positions[p->get_id()] = p->center_of_mass();
			orientations[p->get_id()] = p->orientation();
		}

		comm_buff.send_and_receive_no_wait(true);
		comm_buff.finalize_send_and_receive();

		char * it = comm_buff.recv_buffer_begin();
		while(it != comm_buff.recv_buffer_end()) {
			plint id;
			utils::unpack(it, id);
			utils::unpack(it, positions[id]);
			utils::unpack(it, orientations[id]);
		}

		// Add particles to the immersed boundary object
		for(plint i = 0; i < num_rbcs; ++i) {
			RBCParticle<T> rbc2 = rbc;
			Transform<T> transform;
			transform.translate(-rbc2.center_of_mass())
					 .rotate(orientations[i])
					 .translate(positions[i]);
			rbc2.transform_vertices(transform);
			fsi.add_particle(&rbc2);
		}

		for(plint i = num_rbcs; i < num_particles; ++i) {
			SemiRigidParticle3D<T> platelet2 = platelet;
			Transform<T> transform;
			transform.translate(-platelet2.center_of_mass())
					 .rotate(orientations[i])
					 .translate(positions[i]);
			platelet2.transform_vertices(transform);
			fsi.add_particle(&platelet2);
		}
	}
}

Time create_time()
{
	Time it(0);
	if(load_checkpoint) {
		it.read_checkpoint("checkpoint_time");
	}
	return it;
}

void write_checkpoint(
		const Time & it,
		const ImmersedBoundaryDynamics3D<double, DESCRIPTOR, Periodicity> & fsi,
		MultiBlockLattice3D<double,DESCRIPTOR> & lattice
)
{
	pcout << "Writing checkpoint at iteration " << it << std::endl;
	parallelIO::save(lattice, "checkpoint_lbm", false);
	fsi.save_checkpoint("checkpoint_fsi");
	it.write_checkpoint("checkpoint_time");
}


int main(int argc, char ** argv)
{
	plbInit(&argc, &argv);

	// Read input
	std::string parameter_file;
	try {
		global::argv(1).read(parameter_file);
	} catch(PlbIOException& exception) {
		// Print the corresponding error message.
		pcout << exception.what() << std::endl;
		return EXIT_FAILURE;
	}

	// Read parameters
	ParticleShapeLibrary<double> shape_library;
	read_parameters(shape_library, parameter_file);

	global::IOpolicy().activateParallelIO(true);
	global::directories().setOutputDir(output_folder);
	global::directories().setInputDir(checkpoint_folder);
	io::mkdir(global::directories().getOutputDir().c_str());

	// Data containers
	MultiBlockLattice3D<double,DESCRIPTOR> * lattice = 0;
	GuoRigidWallBoundary<double, DESCRIPTOR> * boundary_condition = 0;

	// Setup lattice and boundary conditions
	setup_flow<double>(lattice, boundary_condition);
	Box3D domain = lattice->getBoundingBox();

	// Initialize immersed boundary module
	ImmersedBoundaryDynamics3D<double, DESCRIPTOR, Periodicity> fsi(
				*lattice,
				shape_library);

	// Add wall interaction
	if(wall_interaction)
		fsi.add_force_decorator(wall_interaction);

	// Add particle-particle interation
	if(pp_interaction)
		fsi.add_pp_force(pp_interaction);

	// Create particles
	create_particles<double>(fsi, *lattice, shape_library);
	fsi.init();

	// Velocity and concentration fields
	std::auto_ptr<MultiTensorField3D<double, 3> > velocity = generateMultiTensorField<double, 3>(domain, extendedEnvelopeWidth);
	std::auto_ptr<MultiScalarField3D<double> > H = generateMultiScalarField<double>(domain, extendedEnvelopeWidth);
	std::auto_ptr<MultiTensorField3D<double, 3> > forceField = generateMultiTensorField<double, 3>(domain, extendedEnvelopeWidth);

	pcout << "===============" << std::endl;
	pcout << "Starting coupled ib-lbm iterations" << std::endl;

	Time it = create_time();
	bool first_iteration = true;

	// Main computation loop
	for(; it <= Nt; ++it) {
		if(it.is_multiple_of(100)) {
			pcout << "Iteration " << it << " / " << Nt << " total volume: " << fsi.compute_total_volume() << std::endl;
		}

		if(it.is_multiple_of(checkpoint_interval) /*&& !first_iteration*/)
			write_checkpoint(it, fsi, *lattice);

		computeVelocity(*lattice, *velocity, domain);

		if(it.is_multiple_of(particle_lw_interval)) {
			pcout << "  Writing lightweight particle data at iteration " << it << std::endl;
			fsi.write_lightweight_particle_data(it.to_plint());
		}

		if(it.is_multiple_of(particle_vtk_interval)) {
			pcout << "  Writing particle vtks at iteration " << it << std::endl;
			fsi.write_particles_as_vtk(it.to_plint());
		}

		if(it.is_multiple_of(flow_vtk_interval)) {
			pcout << "  Writing flow field at iteration " << it << std::endl;

			// Get concentration and force field
			applyProcessingFunctional(new ComputeConcentrationFunctional<double, DESCRIPTOR>, domain, *lattice, *H);
			applyProcessingFunctional(new ForceFieldFunctional<double, DESCRIPTOR>, domain, *lattice, *forceField);

			VtkImageOutput3D<double> vtkOut(createFileName("vtk", it.to_plint(), 6), 1.);
			vtkOut.writeData<3,float>(*velocity, "velocity", 1.);
			vtkOut.writeData<float>(*H, "concentration", 1.);
			vtkOut.writeData<3, float>(*forceField, "force", 1.);
		}

		// Perform IB step
		setExternalVector(*lattice, domain, DESCRIPTOR<double>::ExternalField::forceBeginsAt, external_force);
		applyProcessingFunctional(wrap_ibm_dynamics3D(fsi), domain, *lattice, *velocity);

		// Do lattice boltzmann step
		Profile::start_timer("lbm");
		lattice->collideAndStream();
		Profile::stop_timer("lbm");

		// Just a flag
		first_iteration = false;
	}

	Profile::write_report("profile");

	//parallelIO::save(*lattice, "checkpoint_lbm", false);
	//fsi.save_checkpoint("checkpoint_fsi");
	//it.write_checkpoint("checkpoint_time");

	// Free memory
	delete lattice;
	if(boundary_condition)
		delete boundary_condition;
	if(wall_interaction)
		delete wall_interaction;
	if(pp_interaction)
		delete pp_interaction;
}

