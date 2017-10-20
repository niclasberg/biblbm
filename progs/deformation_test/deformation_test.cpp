/*
 * deformation_test.cpp
 *
 *  Created on: Jun 2, 2015
 *      Author: niber
 */

#include "palabos3D.h"
#include "palabos3D.hh"
#include "fsi/headers.h"
#include "fsi/headers.hh"
#include "../common/calibrate_rbc_shape.h"
#include "../common/rbc_parameters.h"

using namespace plb;
using namespace fsi;

#define DESCRIPTOR descriptors::ForcedPhaseD3Q19Descriptor

// Define the periodicity of the domain (periodicity_x, periodicity_y, periodicity_z)
typedef Periodicity3D<double, true, true, true> Periodicity;

template<class T>
class OpticalTweezerForce : public ForceDecorator<T> {
public:
	void add_id(plint id) { ids_.push_back(id); }
	plint num_ids() const { return ids_.size(); }

	void set_force(const Array<T, 3> & force) { force_ = force; }

	virtual void apply_force(ParticleBase3D<T> * particle) {
		for(std::vector<plint>::iterator it = ids_.begin(); it != ids_.end(); ++it)
			particle->get_node(*it).force += force_;
	}

private:
	Array<T, 3> force_;
	std::vector<plint> ids_;
};

template<class T>
struct XIsLess {
	XIsLess(const RBCParticle<T> & rbc) : rbc_(rbc) { }
	bool operator()(plint i0, plint i1) {
		return rbc_.get_node(i0).pos[0] < rbc_.get_node(i1).pos[0];
	}
private:
	const RBCParticle<T> & rbc_;
};

int main(int argc, char ** argv)
{
	plbInit(&argc, &argv);
	global::IOpolicy().activateParallelIO(true);

	// Read parameter file name	
	std::string parameter_file;
	try {
		global::argv(1).read(parameter_file);
	} catch(PlbIOException& exception) {
		pcout << exception.what() << std::endl;
		return EXIT_FAILURE;
	}

	// Parse parameters
	std::string folder;
	double G, K, K_bend, bend_energy, poisson_ratio, lambda;
	double tau, dx, nu_si, rho_si, deflation_factor; 
	double min_force_si, max_force_si, node_fraction;
	plint Nrelax, force_steps;
	std::string mesh_file;
	try {
		// Open the XMl file
		XMLreader xmlFile(parameter_file.c_str());

		// Folders
		xmlFile["output"]["folder"].read(folder);
		if(folder.at(folder.size() -1) != '/')
			folder.push_back('/');

		// Membrane properties
		xmlFile["rbc"]["shear_modulus"].read(G);
		xmlFile["rbc"]["poisson_ratio"].read(poisson_ratio);
		K = (poisson_ratio + 1)/(1 - poisson_ratio) * G;
		xmlFile["rbc"]["bending_energy"].read(bend_energy);
		K_bend = 2. / std::sqrt(3.) * bend_energy;
		xmlFile["rbc"]["lambda"].read(lambda);

		// RBC forcing
		xmlFile["forcing"]["min_force"].read(min_force_si);
		xmlFile["forcing"]["max_force"].read(max_force_si);
		xmlFile["forcing"]["node_fraction"].read(node_fraction);
		xmlFile["forcing"]["Nrelax"].read(Nrelax);
		xmlFile["forcing"]["force_steps"].read(force_steps);
		
		// Mesh file
		xmlFile["rbc"]["mesh"].read(mesh_file);
		xmlFile["rbc"]["deflation_factor"].read(deflation_factor);

		// Numerics
		xmlFile["units"]["tau"].read(tau);
		xmlFile["units"]["nu_plasma"].read(nu_si);
		xmlFile["units"]["rho_plasma"].read(rho_si);
		xmlFile["units"]["dx"].read(dx);

	} catch (PlbIOException& exception) {
		pcout << exception.what() << std::endl;
		return false;
	}

	global::directories().setOutputDir(folder);
	global::directories().setInputDir(folder);
	io::mkdir(global::directories().getOutputDir().c_str());

	// Load mesh
	ParticleShapeLibrary<double> shape_library;
	//shape_library.read_and_store_mesh("sphere_TRImesh_626n_1248el.msh", "RBC");
	shape_library.read_and_store_mesh(mesh_file.c_str(), "RBC");
	ParticleShape<double> * shape = shape_library.get_by_tag("RBC");

	// Setup unit converter
	double nu_lb = (tau - 0.5) / 3.0;

	UnitConverter<double> unit_converter;
	unit_converter.set_si_to_lb_density_ratio(rho_si/1.);
	unit_converter.set_si_to_lb_length_ratio(dx);
	unit_converter.set_si_to_lb_time_ratio(dx*dx*nu_lb/nu_si);

	RBCParameters<double> rbc_params = create_rbc_params(shape, 
		unit_converter.si_force_per_length_to_lb(G),
		unit_converter.si_force_per_length_to_lb(K),
		unit_converter.si_energy_to_lb(K_bend));

	// Create RBC
	RBCParticle<double> rbc(shape_library.get_by_tag("RBC"), rbc_params);

	// Generate equilibrium shape
	shrink_rbc_volume(rbc, (double) deflation_factor*rbc.shape()->get_volume(), 100000);
	rbc.set_center_of_mass(Array<double, 3>(0, 0, 0));
	rbc.set_minor_axis_orientation(Array<double, 3>(0, 0, 1));

	// Evaluate fluid box dimensions
	double diameter = std::max(rbc.bounding_box().x1 - rbc.bounding_box().x0,
							 std::max(rbc.bounding_box().y1 - rbc.bounding_box().y0,
									  rbc.bounding_box().z1 - rbc.bounding_box().z0));
	int N = std::ceil(diameter);
	Box3D dom(-2*N, 2*N, -N, N, -N, N);

	// LBM relaxation time
	double tau_inner = 0.5 + lambda*(tau - 0.5);
	double omega = 1. / tau;
	double omega_inner = 1. / tau_inner;

	// Print properties
	pcout << "RBC Parameters:" << std::endl;
	pcout << "  Shear modulus: " << rbc_params.G() << std::endl;
	pcout << "  Area modulus " << rbc_params.K() << std::endl;
	pcout << "  Young's modulus: " << rbc_params.youngs_modulus() << std::endl;
	pcout << "  Poisson ratio " << rbc_params.poisson_ratio() << std::endl;
	pcout << "  Relative bending stiffness (kb/GA): " << rbc_params.k_bend / (rbc_params.G() * rbc.area()) << std::endl;
	pcout << "==============" << std::endl;
	pcout << "  Global area stiffness: " << rbc_params.k_area_global << std::endl;
	pcout << "  Local area stiffness: " << rbc_params.k_area_local << std::endl;
	pcout << "  Global volume stiffness: " << rbc_params.k_volume << std::endl;
	pcout << "  Bending stiffness: " << rbc_params.k_bend << std::endl;
	pcout << "  Equilibrium angle: " << rbc_params.theta0 << std::endl;
	pcout << "==============" << std::endl;
	pcout << "Fluid properties:" << std::endl;
	pcout << "  Domain x:" << dom.x0 << ", " << dom.x1 << std::endl;
	pcout << "         y:" << dom.y0 << ", " << dom.y1 << std::endl;
	pcout << "         z:" << dom.z0 << ", " << dom.z1 << std::endl;
	pcout << "  tau (outer): " << tau << std::endl;
	pcout << "  tau (inner): " << tau_inner << std::endl;
	pcout << "==============" << std::endl;
	pcout << "Scales:" << std::endl;
	pcout << "  LB delta x: " << unit_converter.lb_length_to_si(1.) << std::endl;
	pcout << "  LB delta t: " << unit_converter.lb_time_to_si(1.) << std::endl;
	pcout << "==============" << std::endl;

	// Create fluid lattice
	MultiBlockLattice3D<double, DESCRIPTOR> * lattice = generateMultiBlockLattice<double, DESCRIPTOR> (
			dom,
			new GuoExternalForceVOFDynamics<double, DESCRIPTOR>(omega, omega_inner),
			2).release();

	Box3D domain = lattice->getBoundingBox();

	lattice->toggleInternalStatistics(false);
	lattice->periodicity().toggle(0, Periodicity::get_x());
	lattice->periodicity().toggle(1, Periodicity::get_y());
	lattice->periodicity().toggle(2, Periodicity::get_z());

	setExternalVector(*lattice, domain, DESCRIPTOR<double>::ExternalField::forceBeginsAt, Array<double, 3>(0., 0., 0.));

	initializeAtEquilibrium(*lattice, domain, 1., Array<double, 3>(0, 0, 0));

	// Create immersed boundary
	ImmersedBoundaryDynamics3D<double, DESCRIPTOR, Periodicity> fsi(
			*lattice,
			shape_library );
	std::auto_ptr<MultiTensorField3D<double, 3> > velocity = generateMultiTensorField<double, 3>(domain, 2);

	// Find vertices to apply the force to
	OpticalTweezerForce<double> negativeXForce, positiveXForce;
	plint numNodes = std::ceil(node_fraction * rbc.count_nodes());

	{
		std::vector<plint> nodeIds;
		for(plint i = 0; i < rbc.count_nodes(); ++i)
			nodeIds.push_back(i);

		// Sort the ids according to x-position in the RBC particle
		std::sort(nodeIds.begin(), nodeIds.end(), XIsLess<double>(rbc));
		for(plint i = 0; i < numNodes; ++i) {
			negativeXForce.add_id(nodeIds[i]);
			positiveXForce.add_id(nodeIds[nodeIds.size() - i - 1]);
		}
	}

	std::cout << "Marked " << negativeXForce.num_ids() << " vertices for negative x force and " << positiveXForce.num_ids() << " for positive" << std::endl;

	fsi.add_force_decorator(&negativeXForce);
	fsi.add_force_decorator(&positiveXForce);

	fsi.add_particle(&rbc);
	fsi.init();

	double minForce = unit_converter.si_force_to_lb(min_force_si);
	double maxForce = unit_converter.si_force_to_lb(max_force_si);

	computeVelocity(*lattice, *velocity, domain);

	// Main computation loop
	std::stringstream ss;
	ss << global::directories().getOutputDir() << "deformation" << global::mpi().getRank() << ".txt";
	std::fstream fout(ss.str().c_str(), std::ios::out);
	if(!fout.good())
		std::cerr << "Could not open file for writing" << std::endl;

	int it = 0;

	for(int fIt = 1; fIt <= force_steps; ++fIt) {
		double totalForce = minForce + (maxForce - minForce) * (double) fIt / (double) (force_steps);
		negativeXForce.set_force(Array<double, 3>(-totalForce / (double) negativeXForce.num_ids(), 0., 0.));
		positiveXForce.set_force(Array<double, 3>(totalForce / (double) positiveXForce.num_ids(), 0., 0.));

		pcout << "Force iteration " << fIt << ", total force = " << totalForce << std::endl;

		{
			VtkImageOutput3D<double> vtkOut(createFileName("vtk", fIt, 4), 1.);
			vtkOut.writeData<3,float>(*velocity, "velocity", 1.);
			fsi.write_particles_as_vtk(fIt);
		}

		for(int subIt = 0; subIt <= Nrelax; ++subIt) {
			if(subIt == Nrelax || (subIt % 1000) == 0) {
				// Compute diameters
				ParticleBase3D<double> * p;
				double Dx, Dyz;
				Dx = 0.;
				Dyz = 0.;
				if(fsi.get_particle(0, p)) {
					for(plint i = 0; i < p->count_nodes(); ++i) {
						Array<double, 3> d = p->get_node(i).pos - p->center_of_mass();
						Dx = std::max(Dx, 2*std::abs(d[0]));
						Dyz = std::max(Dyz, 2*std::sqrt(d[1]*d[1] + d[2]*d[2]));
					}
					fout << it << ", " << unit_converter.lb_force_to_si(totalForce) << ", " <<
						unit_converter.lb_length_to_si(Dx) << ", " <<
						unit_converter.lb_length_to_si(Dyz) << std::endl;

					std::cout << "  Sub iteration " << subIt << " / " << Nrelax << " Dx = " << Dx << ", Dyz = " << Dyz <<  std::endl;
				}
			}

			// Perform IB step
			computeVelocity(*lattice, *velocity, domain);

			// Do lattice boltzmann step
			Profile::start_timer("lbm");
			lattice->collideAndStream();
			Profile::stop_timer("lbm");

			// Do IB step
			setExternalVector(*lattice, domain, DESCRIPTOR<double>::ExternalField::forceBeginsAt, Array<double, 3>(0, 0, 0));
			applyProcessingFunctional(wrap_ibm_dynamics3D(fsi), domain, *lattice, *velocity);
			++it;
		}
	}
	fout.close();
}
