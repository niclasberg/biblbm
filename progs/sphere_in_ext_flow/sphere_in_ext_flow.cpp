// Comment away this line to disable profiling:
//#define FSI_PROFILE

#include "palabos3D.h"
#include "palabos3D.hh"
#include "fsi/headers.h"
#include "fsi/headers.hh"
#include "../common/shear_flow_setup.h"
#include "../common/extensional_flow_setup.h"
#include "../common/particle_position_regulator.h"

namespace plb {
namespace fsi {

// Function object to execute at each timestep before the LBM step
template<class T, template<typename U> class Descriptor, class Periodicity>
class MultiDirectForcingFunctional3D : public BoxProcessingFunctional3D_LT<T, Descriptor, T, 3> {
public:
	MultiDirectForcingFunctional3D(ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity> & ibm, plint relaxIterations)
	: ibm_(ibm), relaxIterations_(relaxIterations)
	{ 

	}

	virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice, TensorField3D<T,3> & velocity) {
		if( ! ibm_.is_init()) {
			pcerr << "The immersed boundary module has not been initialized! Call the ImmersedBoundaryDynamics3D<...>::init() " <<
					"method after all particles have been added!" << std::endl;
			exit(-1);
		}

		Dot3D offset = computeRelativeDisplacement(lattice, velocity);
		
		// Set the force on the fluid phase to zero
		for(plint i = domain.x0; i <= domain.x1; ++i)
			for(plint j = domain.y0; j <= domain.y1; ++j)
				for(plint k = domain.z0; k <= domain.z1; ++k) {
					T * g = lattice.get(i, j, k).getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
					g[0] = 0.;
					g[1] = 0.;
					g[2] = 0.;
				}

		// Container to store the total node forces
		std::vector<Array<T, 3> > forces;

		// Check if the current processor owns the particle
		bool hasParticle = false;
		ParticleBase3D<T> * p;
		if(ibm_.get_particle(0, p)) {
			hasParticle = true;
			forces.resize(p->count_nodes());
			for(plint i = 0; i < forces.size(); ++i)
				forces[i].resetToZero();
		}
 
		for(plint it = 0; it < this->relaxIterations_; ++it) {
			// Compute velocity
			// In the Guo formulation this is given by
			//  rho u = sum_i f_i e_i + F/2
			// where F is the sum of external forces
			for(plint i = domain.x0; i <= domain.x1; ++i)
				for(plint j = domain.y0; j <= domain.y1; ++j)
					for(plint k = domain.z0; k <= domain.z1; ++k)
						lattice.get(i, j, k)
							.computeVelocity(velocity.get(i+offset.x, j+offset.y, k+offset.z));

			// Interpolate the velocity to the particles
			ibm_.interpolate_velocity(domain, velocity);

			// Compute fsi force and spread to the fluid nodes
			ibm_.compute_and_spread_forces(domain, lattice);

			if(hasParticle) {
				Array<T, 3> totalForce; totalForce.resetToZero();
				for(plint j = 0; j < p->count_nodes(); ++j) {
					forces[j] += p->get_node(j).force;
					totalForce += p->get_node(j).force;
				}
				//std::cout << " Subiteration " << it << ", total force norm = " << norm(totalForce) << std::endl;
			}
		}

		// Set the particle node forces to the total
		if(hasParticle) {
			for(plint j = 0; j < p->count_nodes(); ++j) {
				p->get_node(j).force = -forces[j];
			}
		}

		// Integrate the equations of motion for the particle
		ibm_.move_vertices();
	}

	virtual MultiDirectForcingFunctional3D * clone() const
	{
		return new MultiDirectForcingFunctional3D(*this);
	}

	virtual BlockDomain::DomainT appliesTo() const {
		return BlockDomain::bulk;
	}

	virtual void getModificationPattern(std::vector<bool>& isWritten) const {
		isWritten[0] = true;
		isWritten[1] = true;
	}

	virtual void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
		modified[0] = modif::staticVariables;
		modified[1] = modif::staticVariables;
	}
private:
	ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity> & ibm_;
	plint relaxIterations_;
};
}
}

using namespace plb;
using namespace fsi;
#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor

// Define the periodicity of the domain (periodicity_x, periodicity_y, periodicity_z)
typedef Periodicity3D<double, true, false, true> Periodicity;

int main(int argc, char ** argv)
{
	// Initialize palabos
	plbInit(&argc, &argv);

	// Setup input/output
	global::IOpolicy().activateParallelIO(true);

	// Numerics
	double relaxationFactor = 0.3;	// relaxation factor for the forcing iterations
	plint relaxIterations = 4;

	// Domain dimensions
	plint geoNx = 96, geoNy = 96, geoNz = 96;

	// Number of iterations
	plint Nt = 1200001;

	// Number of iterations between each write
	plint output_interval = 200;
	plint particle_data_output_interval = 10;

	// Input parameters
	double tau = 0.55; 			// LBM relaxation time
	double Re_p = 1;			// Particle Reynolds number (= (radius)^2 * gamma / nu)
	double rho_p = 1.;			// Particle density
	const char * mesh_file = "sphere_r8.gmsh";
	double nu_lb = (tau - 0.5) / 3.0;
	double omega = 1. / tau;

	// Create output folder
	std::stringstream ss;
	ss << "sphereData" << geoNx << "x" << geoNy << "x" << geoNz << "Re" << Re_p << "tau" << tau << "/";
	std::string folder = ss.str();
	io::mkdir(folder.c_str());
	global::directories().setOutputDir(folder);
	global::directories().setInputDir(global::directories().getOutputDir());

	// Read and store the mesh
	ParticleShapeLibrary<double> shape_library;
	shape_library.read_and_store_mesh(mesh_file, "sphere");
	ParticleShape<double> * shape = shape_library.get_by_tag("sphere");

	// Deduce flow parameters
	double a_lb = std::sqrt(shape->get_area() / (4. * M_PI));
	double gamma_lb = Re_p * nu_lb / (a_lb*a_lb);
	double uMax = gamma_lb * std::sqrt(geoNy*geoNy + geoNz*geoNz) / 2.0;

	// Create the particle and set its position
	Array<double, 3> particlePosition(0.5*geoNx, 0.8*geoNy, 0.5*geoNz);
	RigidParticle3D<double> particle(shape);
	particle.relaxation_factor() = relaxationFactor;
	particle.density() = rho_p;
	//particle.angular_velocity()[2] = -gamma_lb / 2.0;
	particle.set_center_of_mass(particlePosition);
	particle.velocity() = extensionalFlowVelocity(particlePosition, geoNx, geoNy, gamma_lb);

	pcout << "==============" << std::endl;
	pcout << "Non-dimensional parameters:" << std::endl;
	pcout << "  Re_p: " << gamma_lb * util::sqr(a_lb) / nu_lb << std::endl;
	pcout << "  Re_p: " << gamma_lb * util::sqr(a_lb) / nu_lb << std::endl;
	pcout << "==============" << std::endl;
	pcout << "Numerics:" << std::endl;
	pcout << "  uMax:" << uMax << std::endl;
	pcout << "  tau:" << tau << std::endl;
	pcout << "  forcing relaxation factor: " << particle.relaxation_factor() << std::endl;
	pcout << "  forcing relaxation iterations: " << relaxIterations << std::endl;

	// Create lattice
	MultiBlockLattice3D<double, DESCRIPTOR> lattice(
			geoNx+1, geoNy+1, geoNz+1,
			new GuoExternalForceBGKdynamics<double, DESCRIPTOR>(omega) );
	Box3D domain = lattice.getBoundingBox();
	lattice.toggleInternalStatistics(false);

	OnLatticeBoundaryCondition3D<double,DESCRIPTOR> * boundaryCondition
			= createInterpBoundaryCondition3D<double,DESCRIPTOR>();

	// Setup boundary conditions
	//shearFlowSetup(lattice, uMax, *boundaryCondition);
	extensionalFlowSetup(lattice, gamma_lb, *boundaryCondition);

	// Immersed boundary module
	ImmersedBoundaryDynamics3D<double, DESCRIPTOR, Periodicity> fsi(
			lattice,
			shape_library);

	// Velocity and concentration fields
	std::auto_ptr<MultiTensorField3D<double, 3> > velocity = generateMultiTensorField<double, 3>(domain);

	// Add particle
	fsi.add_particle(&particle);

	// Initialize fsi
	fsi.init();

	// Open file for particle data output for the current processor
	std::string outputFileName;
	{
		std::stringstream ss;
		ss << folder << "particleData" << global::mpi().getRank() << ".csv";
		outputFileName = ss.str();
	}

	// Check if file exists
	bool fileExisted = std::ifstream(outputFileName.c_str()).good();

	std::ofstream lightweightOut(outputFileName.c_str(), std::ios_base::app);
	if( ! fileExisted) {
		// If file did not already exist, write header
		lightweightOut << "iteration,x,y,z,ux,uy,uz,Fx,Fy,Fz,Tx,Ty,Tz" << std::endl;
	}

	pcout << "===============" << std::endl;
	pcout << "Starting coupled ib-lbm iterations" << std::endl;
	for(plint it = 0; it <= Nt; ++it) {
		if((it % 100) == 0) {
			ParticleBase3D<double> * p;
			if(fsi.get_particle(0, p)) {
				RigidParticle3D<double> * pp = dynamic_cast<RigidParticle3D<double> *>(p);
				std::cout << "Iteration " << it << " / " << Nt 
						<< ", particle velocity = (" << pp->velocity()[0] << ", " << pp->velocity()[1] << ", " << pp->velocity()[2]
						<< "), angular velocity = (" << pp->angular_velocity()[0] << ", " << pp->angular_velocity()[1] << ", " << pp->angular_velocity()[2] << ")" << std::endl;
				//std::cout << "omega / gamma = " << (- pp->angular_velocity()[2] / gamma_lb) << std::endl;
			}
		}
		
		// Write lightweight particle data
		if((it % particle_data_output_interval) == 0) {
			ParticleBase3D<double> * pBase;
			if(fsi.get_particle(0, pBase)) {
				// Cast to rigid particle
				RigidParticle3D<double> * p = dynamic_cast<RigidParticle3D<double> *>(pBase);
				if( ! p) {
					// This should not happen since the particle should be of type RigidParticle3D
					std::cout << "Error while while casting particle type to RigidParticle3D" << std::endl;
					exit(-1);
				}

				// Write data
				lightweightOut << it << ","
					<< p->center_of_mass()[0] << "," << p->center_of_mass()[1] << "," << p->center_of_mass()[2] << "," 
					<< p->velocity()[0] << "," << p->velocity()[1] << "," << p->velocity()[2] << ","
					<< p->force()[0] << "," << p->force()[1] << "," << p->force()[2] << ","
					<< p->torque()[0] << "," << p->torque()[1] << "," << p->torque()[2] << std::endl;
			}
		}

		/*if((it % output_interval) == 0) {
			fsi.write_lightweight_particle_data(it);
		}*/

		if((it % output_interval) == 0) {
			computeVelocity(lattice, *velocity, domain);
			fsi.write_particles_as_vtk(it);

			VtkImageOutput3D<double> vtkOut(createFileName("vtk", it, 6), 1.);
			vtkOut.writeData<3,float>(*velocity, "velocity", 1.);
		}

		// Perform FSI stuff
		applyProcessingFunctional(
			new MultiDirectForcingFunctional3D<double, DESCRIPTOR, Periodicity>(fsi, relaxIterations), 
			domain, 
			lattice, *velocity);

		// Do LBM
		lattice.collideAndStream();
	}

	//Profile::write_report("profile");

	parallelIO::save(lattice, "checkpoint_lbm", false);
	fsi.save_checkpoint("checkpoint_fsi");

	global::mpi().barrier();
}

