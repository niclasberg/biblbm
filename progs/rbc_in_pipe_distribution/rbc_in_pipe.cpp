//#define FSI_PROFILE

#include "palabos3D.h"
#include "palabos3D.hh"
#include "fsi/headers.h"
#include "fsi/headers.hh"
#include "../common/calibrate_rbc_shape.h"
#include "../common/util.h"
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
plint margin = 2;
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
Array<double, 3> external_force; 	// Pressure gradient
PipeBoundary<double> pipe_boundary;
bool second_order_bcs = true;		// true  => use Guo extrapolation boundary conditions for cylinder surface
									// false => use staircase bounce back at cylinder surface
std::string output_folder;
std::string checkpoint_folder;
PPInteractionForce<double> * pp_interaction = 0;
ForceDecorator<double> * wall_interaction = 0;

struct ParticleTmp {
	ParticleBase3D<double> * particle;
	std::string name;
	double numParticles;
};
std::vector<ParticleTmp> particleTypes;


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

bool read_parameters(ParticleShapeLibrary<double> & shape_library, std::string file_name)
{
	// Input parameters
	double tau; 			// Relaxation parameter
	double Re_d;			// Channel Reynolds number
	double lambda;

	try {
		// Open the XMl file
		XMLreader xmlFile(file_name.c_str());

		// Pipe geometry
		pcout << "Reading geometry" << std::endl;
		xmlFile["geometry"]["diameter"].read(diameter);
		xmlFile["geometry"]["length"].read(length);

		// Set boundary dimensions
		pipe_boundary.y0 = 0.;
		pipe_boundary.z0 = 0.;
		pipe_boundary.radius = (double) diameter / 2.;
		pipe_boundary.length = (double) length;
		double volume = M_PI * util::sqr(pipe_boundary.radius) * pipe_boundary.length;

		// Read physical parameters
		pcout << "Reading physical parameters" << std::endl;
		xmlFile["physics"]["Re_d"].read(Re_d);
		xmlFile["physics"]["lambda"].read(lambda);

		// Read numerics
		pcout << "Reading numberical parameters" << std::endl;
		xmlFile["numerics"]["tau"].read(tau);
		xmlFile["numerics"]["max_iterations"].read(Nt);

		// Compute fluid viscosity and relaxation frequencies
		double nu_lb = (tau - 0.5) / 3.0;
		double tau_inner = 0.5 + lambda*(tau - 0.5);
		omega = 1. / tau;
		omega_inner = 1. / tau_inner;

		// Compute maximal velocity ( = 2 * u_mean) and effective shear rate
		uMax = 2. * nu_lb * Re_d / (double) diameter;
		double gamma_lb = 4. * uMax / (double) diameter;

		// External force
		external_force[0] = 4. * nu_lb * gamma_lb / (double) diameter;
		external_force[1] = 0.;
		external_force[2] = 0.;

		pcout << "Numerics:" << std::endl;
		pcout << "  uMax:" << uMax << std::endl;
		pcout << "  tau (outer fluid):" << tau << std::endl;
		pcout << "  tau (inner fluid):" << tau_inner << std::endl;
		pcout << "==============" << std::endl;

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

		// Read particle data
		XMLreaderProxy particleParams = xmlFile["particles"];
		for(int i = 0; i < particleParams.getChildren().size(); ++i) {
			XMLreaderProxy node = particleParams.getChildren()[i];
			std::string type = node.getName();
			std::string subType;
			node["type"].read(subType);
			pcout << type << " - " << subType << std::endl;

			// Wall interaction
			if(type.compare("wall_interaction") == 0) {
				if(subType.compare("morse") == 0) {
					wall_interaction =  new WallInteraction<double, MorsePotential<double> >(
							pipe_boundary,
							create_morse_potential<double>(xmlFile["particles"]["wall_interaction"]));
				} else if(subType.compare("spring") == 0) {
					wall_interaction =  new WallInteraction<double, SpringPotential<double> >(
							pipe_boundary,
							create_spring_potential<double>(xmlFile["particles"]["wall_interaction"]));
				} else if(subType.compare("none") == 0) {

				} else {
					std::cerr << "Unknown wall interaction type " << subType << std::endl;
				}
			} else if(type.compare("particle_interaction") == 0) {
				if(subType.compare("morse") == 0) {
					pp_interaction =  new PPPotentialForce<double, MorsePotential<double> >(
							create_morse_potential<double>(xmlFile["particles"]["particle_interaction"]));
				} else if(subType.compare("spring") == 0) {
					pp_interaction =  new PPPotentialForce<double, SpringPotential<double> >(
							create_spring_potential<double>(xmlFile["particles"]["particle_interaction"]));
				} else if(subType.compare("none") == 0) {

				} else {
					std::cerr << "Unknown particle interaction type " << subType << std::endl;
					return false;
				}
			} else {
				// Assume that it is a particle
				// Read mesh
				std::string meshFile;
				node["mesh"].read(meshFile);
				shape_library.read_and_store_mesh(meshFile.c_str(), type.c_str());
				ParticleShape<double> * shape = shape_library.get_by_tag(type.c_str());

				// If we're starting from a checkpoint, the particles in the domain will already be initialized
				if(!load_checkpoint) {
					double Ca;
					double concentration;

					// Read species concentration
					node["concentration"].read(concentration);

					// Temporary particle object
					ParticleTmp particleTmp;
					particleTmp.name = type;

					if(subType.compare("rbc") == 0) {
						// Read rbc properties
						double k_bend_rel;
						double poisson_ratio;
						double deflationFactor;

						node["Ca"].read(Ca);
						node["poisson_ratio"].read(poisson_ratio);
						node["Xi_b"].read(k_bend_rel);
						node["mesh_deflation"].read(deflationFactor);

						// Create particle properties
						double a_lb = std::sqrt(shape->get_area() / (4. * M_PI));

						// Membrane parameters
						double G_lb = a_lb * gamma_lb * nu_lb / Ca;
						double K_area = (poisson_ratio + 1)/(1 - poisson_ratio) * G_lb;
						double K_bend = k_bend_rel * G_lb * shape->get_area();

						// Create particle initializer
						RBCShapeInitializer<double> rbc_initializer(new RBCParticle<double>(shape));

						// Create equilbrium shape
						rbc_initializer.compute_rbc_parameters(G_lb, K_area, K_bend);
						rbc_initializer.shrink_volume(deflationFactor, 40000);

						// Store particle
						particleTmp.particle = rbc_initializer.get_rbc().clone();

					} else if(subType.compare("platelet") == 0) {
						// Read properties
						double k_bend_rel;

						node["Ca"].read(Ca);
						node["Xi_b"].read(k_bend_rel);

						// Create particle properties
						double a_lb = std::sqrt(shape->get_area() / (4. * M_PI));

						// Membrane parameters
						double G_lb = a_lb * gamma_lb * nu_lb / Ca;
						double K_bend = k_bend_rel * G_lb * shape->get_area();

						// Store particle
						RBCParticle<double> * platelet = new RBCParticle<double>(shape, create_platelet_params(shape, G_lb, 0., K_bend));
						platelet->set_should_voxelize(false);
						particleTmp.particle = platelet;

					} else if(subType.compare("semirigid") == 0) {
						// Read platelet properties
						node["Ca"].read(Ca);

						// Compute average link length
						SemiRigidParticleParams<double> platelet_params;
						double lavg = 0;
						for(ParticleShape<double>::link_const_iterator it = shape->links_begin();
								it != shape->links_end(); ++it)
							lavg += it->length;
						lavg /= shape->count_links();
						platelet_params.l0 = lavg;

						// Compute shear modulus
						double a_p_lb = std::sqrt(shape->get_area() / (4. * M_PI));
						double G_lb = a_p_lb * gamma_lb * nu_lb / Ca;

						// Set in plane energy (proportional to the shear modulus)
						platelet_params.k_in_plane = 1.;
						platelet_params.k_in_plane = G_lb / platelet_params.shear_modulus();
						platelet_params.k_out_of_plane = 10*platelet_params.k_in_plane;

						// Create semirigid particle
						SemiRigidParticle3D<double> * particle = new SemiRigidParticle3D<double>(shape, platelet_params);
						particle->set_center_of_mass(Array<double, 3>(0, 0, 0));
						particle->set_should_voxelize(false);
						particleTmp.particle = particle;

					} else {
						std::cerr << "Unknown particle type " << subType << std::endl;
						return false;
					}
					particleTmp.numParticles = std::ceil(concentration *  volume / particleTmp.particle->volume());
					particleTypes.push_back(particleTmp);
				}
			}
		}
	} catch (PlbIOException& exception) {
		pcout << exception.what() << std::endl;
		return false;
	}

	/*pcout << "RBC Parameters:" << std::endl;
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
	*/

	return true;
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
	if( ! second_order_bcs) {
		defineDynamics<T,DESCRIPTOR>(*fluidLattice,
				fluidLattice->getBoundingBox(),
				new BoundarySurface<T>(pipe_boundary),
				new BounceBack<T,DESCRIPTOR>() );
	}

	// Optimization: turn off dynamics for the cells outside of the cylindrical domain
	/*defineDynamics<T,DESCRIPTOR>(
				*fluidLattice,
				fluidLattice->getBoundingBox(),
				new BoundaryOutside<T>(pipe_boundary, std::sqrt(2)),
				new plb::NoDynamics<T,DESCRIPTOR>() );*/

	// Read checkpoint file if requested
	if(load_checkpoint) {
		parallelIO::load("checkpoint_lbm", *fluidLattice, false);
	} else {
		setExternalVector(*fluidLattice,
				fluidLattice->getBoundingBox(),
				DESCRIPTOR<T>::ExternalField::forceBeginsAt,
				external_force);
		setExternalScalar(*fluidLattice,
				fluidLattice->getBoundingBox(),
				DESCRIPTOR<T>::ExternalField::phaseBeginsAt,
				0.);

		// Default-initialize at equilibrium with zero-velocity and fixed pressure.
		initializeAtEquilibrium(*fluidLattice,
				fluidLattice->getBoundingBox(),
				PipeFlowDensityAndVelocity<double>(pipe_boundary, uMax));

		fluidLattice->initialize();

		pipe_boundary.writeVTK(global::directories().getOutputDir() + std::string("boundary.vtu"));
	}

	if(second_order_bcs) {
		boundary_condition = new GuoRigidWallBoundary<T, DESCRIPTOR>(pipe_boundary, *fluidLattice);
		boundary_condition->insert();

		fluidLattice->executeInternalProcessors();
	}

	if(! load_checkpoint)
		parallelIO::save(*fluidLattice, "checkpoint_lbm", false);
}

// Create an initial particle distribution by initializing the
// particles randomly at a fraction of their size. The particles
// are then grown to their full size while under the influence of
// repulsive forces.
template<class T>
void create_particles(
		ImmersedBoundaryDynamics3D<T, DESCRIPTOR, Periodicity> & fsi,
		MultiBlockLattice3D<T,DESCRIPTOR> & lattice)
{
	// Create intial condition for particles
	if(load_checkpoint) {
		pcout << "Loading particle checkpoint" << std::endl;
		fsi.load_checkpoint("checkpoint_fsi");
	} else {
		if( ! particleTypes.empty()) {
			plint Nit = 100000;
			T initial_scale = 0.2;

			pcout << "Generating initial distribution" << std::endl;

			// Create a new immersed boundary dynamics object for the rigid body simulation
			ParticleShapeLibrary<T> shape_library2;
			ImmersedBoundaryDynamics3D<T, DESCRIPTOR, Periodicity> fsi2(lattice, shape_library2);

			// Create a boundary that is a little bit smaller than the actual one,
			// this ensures that all particles are within the domain
			PipeBoundary<T> boundary2 = pipe_boundary;
			boundary2.radius = pipe_boundary.radius - 1;

			plint num_particles = 0;
			T max_radius = 0;
			std::vector<RigidParticle3D<T> *> rigidParticles;
			std::vector<T> particleInitialRadii;
			for(int i = 0; i < particleTypes.size(); ++i) {
				// Compute the total number of particles
				num_particles += particleTypes[i].numParticles;

				// Store the equilibrium shape
				particleTypes[i].particle->store_shape(shape_library2, particleTypes[i].name.c_str());
				ParticleShape<T> * rigid_shape = shape_library2.get_by_tag(particleTypes[i].name.c_str());

				// Create a rigid particle with the same shape but 90% smaller
				RigidParticle3D<T> * rigid_particle = new RigidParticle3D<T>(rigid_shape);
				rigid_particle->scale() = initial_scale;
				rigid_particle->density() = 1.;
				rigid_particle->damping() = 0.1;
				rigid_particle->update();
				rigidParticles.push_back(rigid_particle);

				// Update max radius
				T initialRadius = initial_scale*rigid_shape->get_radius() + 1.;
				particleInitialRadii.push_back(initialRadius);
				max_radius = std::max(max_radius, initialRadius);
			}

			// Generate initial distribution
			std::vector<Array<T, 3> > positions;
			std::vector<Quaternion<T> > orientations;
			if(global::mpi().isMainProcessor()) {
				// Generate positions
				ParticlePositionInitializer<T> position_initializer(
						lattice.getBoundingBox(),
						&boundary2, max_radius);

				for(int i = 0; i < particleTypes.size(); ++i)
					position_initializer.generate_points(particleTypes[i].numParticles, particleInitialRadii[i], i);

				// Get positions
				for(int i = 0; i < particleTypes.size(); ++i) {
					std::vector<Array<T, 3> > posTmp;
					position_initializer.get_points(i, posTmp);
					positions.insert(positions.end(), posTmp.begin(), posTmp.end());
				}

				// Generate random orientations
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

			pcout << "Initial distribution created with " << num_particles << " particles at " << 100*initial_scale << "% of their full size" << std::endl;

			if(num_particles == 0) {
				pcout << "No particles to add, skipping initialization step" << std::endl;
			} else {

				pcout << "Growing particles to their full size" << std::endl;

				// Create and insert particles to the fsi object
				{
					int N;
					for(int i = 0, N = 0; i < particleTypes.size(); ++i)
						for(plint j = 0; j < particleTypes[i].numParticles; ++j, ++N) {
							RigidParticle3D<T> * p = rigidParticles[i];
							p->scale() = initial_scale;
							p->set_center_of_mass(positions[N]);
							p->orientation() = orientations[N];
							p->update();
							fsi2.add_particle(p);
						}
				}

				// Initialize fsi
				fsi2.init();

				// Potential force
				SpringPotential<T> potential(1, 0.5);

				// Grow the particles while under the influence of collision forces
				T scale;
				for(plint it = 0; it <= Nit; ++it) {
					if((it % 100) == 0)
						pcout << "  Iteration " << it << " / " << Nit << std::endl;
					//if((it % 500) == 0) fsi2.write_particles_as_vtk(it);

					// Rescale all particles 
					T fraction = (T)it / (T) Nit;
					scale = initial_scale + (1. - initial_scale)*std::pow(fraction, 1./3.);

					if(it == Nit)
						scale = 1.;

					for(typename ImmersedBoundaryDynamics3D<T, DESCRIPTOR, Periodicity>::ObjMapIterator iter = fsi2.particles_begin();
							iter != fsi2.particles_end(); ++iter) {
						plb::fsi::RigidParticle3D<T> * p = dynamic_cast<plb::fsi::RigidParticle3D<T> *>(iter->second);
						if(p) {
							p->scale() = scale;
							p->update();
						}
					}

					fsi2.synchronize_particle_states();
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
				{
					int N;
					for(int i = 0, N = 0; i < particleTypes.size(); ++i) {
						for(int j = 0; j < particleTypes[i].numParticles; ++j, ++N) {
							ParticleBase3D<T> * p = particleTypes[i].particle->clone();
							DeformableParticle3D<T> * pd = dynamic_cast<DeformableParticle3D<T> *>(p);
							RigidParticle3D<T> * pr = dynamic_cast<RigidParticle3D<T> *>(p);
							if(pd) {
								Transform<T> transform;
								transform.translate(-pd->center_of_mass())
										 .rotate(orientations[N])
										 .translate(positions[N]);
								pd->transform_vertices(transform);
							} else if(pr) {
								pr->orientation() = orientations[N];
								pr->set_center_of_mass(positions[N]);
							} else {
								std::cerr << "An error occured during particle type casting" << std::endl;
							}
							p->update();

							fsi.add_particle(p);
							delete p;
						}
					}
				}
			}
		}

		fsi.save_checkpoint("checkpoint_fsi");
	}

	for(int i = 0; i < particleTypes.size(); ++i)
		delete particleTypes[i].particle;
	particleTypes.clear();
}

Time create_time()
{
	Time it(0);
	if(load_checkpoint) {
		it.read_checkpoint("checkpoint_time");
	} else {
		it.write_checkpoint("checkpoint_time");
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
	if( ! read_parameters(shape_library, parameter_file)) {
		std::cerr << "Could not read the input XML file" << std::endl;
		exit(-1);
	}

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

	fsi.set_boundary(&pipe_boundary);

	// Add wall interaction
	if(wall_interaction)
		fsi.add_force_decorator(wall_interaction);

	// Add particle-particle interation
	if(pp_interaction)
		fsi.add_pp_force(pp_interaction);

	// Create particles
	create_particles<double>(fsi, *lattice);
	fsi.init();

	// Velocity and concentration fields
	std::auto_ptr<MultiTensorField3D<double, 3> > velocity = generateMultiTensorField<double, 3>(domain, extendedEnvelopeWidth);
	std::auto_ptr<MultiScalarField3D<double> > H = generateMultiScalarField<double>(domain, extendedEnvelopeWidth);
	std::auto_ptr<MultiTensorField3D<double, 3> > forceField = generateMultiTensorField<double, 3>(domain, extendedEnvelopeWidth);

	// Print out the lattice structure
	{
		std::auto_ptr<MultiScalarField3D<double> > ind = generateMultiScalarField<double>(domain, extendedEnvelopeWidth);
		applyProcessingFunctional(new GuoRigidWallBoundaryDebugger<double, DESCRIPTOR>(*boundary_condition), domain, *ind);
		VtkImageOutput3D<double> vtkOut(std::string("latticeStructure"), 1.);
		vtkOut.writeData<float>(*ind, "boundaryNodes", 1.);
	}

	pcout << "===============" << std::endl;
	pcout << "Starting coupled ib-lbm iterations" << std::endl;

	Time it = create_time();
	bool first_iteration = true;

	// Main computation loop
	for(; it <= Nt; ++it) {
		if(it.is_multiple_of(100)) {
			pcout << "Iteration " << it << " / " << Nt
					<< ", particle volume fraction: "
					<< fsi.compute_total_volume() / (M_PI * util::sqr(pipe_boundary.radius) * pipe_boundary.length) << std::endl;
		}

		if(it.is_multiple_of(checkpoint_interval) && !first_iteration)
			write_checkpoint(it, fsi, *lattice);

		// Compute velocity
		//applyProcessingFunctional(new VelocityComputer3D<double, DESCRIPTOR>(), domain, *lattice, *velocity);
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

		// Do Lattice Boltzmann step
		Profile::start_timer("lbm");
		lattice->collideAndStream();
		Profile::stop_timer("lbm");

		// Just a flag
		first_iteration = false;
		//break;
	}

	Profile::write_report("profile");

	parallelIO::save(*lattice, "checkpoint_lbm", false);
	fsi.save_checkpoint("checkpoint_fsi");
	it.write_checkpoint("checkpoint_time");

	// Free memory
	delete lattice;
	if(boundary_condition)
		delete boundary_condition;
	if(wall_interaction)
		delete wall_interaction;
	if(pp_interaction)
		delete pp_interaction;
	for(plint i = 0; i < particleTypes.size(); ++i)
		delete particleTypes[i].particle;
}

