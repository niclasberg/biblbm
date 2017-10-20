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

int main(int argc, char ** argv)
{
	plbInit(&argc, &argv);

	global::IOpolicy().activateParallelIO(true);
	global::directories().setOutputDir("pow_model_fixed2/");
	global::directories().setInputDir("pow_model_fixed2/");
	io::mkdir(global::directories().getOutputDir().c_str());

	// Load mesh
	ParticleShapeLibrary<double> shape_library;

	//shape_library.read_and_store_mesh("sphere_TRImesh_1250n_2496el.msh", "RBC");
	shape_library.read_and_store_mesh("sphere_TRImesh_626n_1248el.msh", "RBC");
	//shape_library.read_and_store_mesh("sphere_TRImesh_2496n_4988el.msh", "RBC");
	ParticleShape<double> * shape = shape_library.get_by_tag("RBC");

	// Viscosity contrast
	double lambda = 5;

	// Membrane properties
	double G = 5.5e-6;
	double K = 300 * 1.e-3;
	double K_bend_rel;

	global::argv(1).read(K_bend_rel);

	// Setup unit converter
	double dx = 0.5e-6;
	//double a_rbc_si = 4e-6, a_rbc_lb = shape->get_radius();
	double tau = 1;
	double nu_lb = (tau - 0.5) / 3.0, nu_si = 1.46e-6;

	UnitConverter<double> unit_converter;
	unit_converter.set_si_to_lb_density_ratio(1025/1);
	unit_converter.set_si_to_lb_length_ratio(dx);
	unit_converter.set_si_to_lb_time_ratio(dx*dx*nu_lb/nu_si);

	double G_lb = unit_converter.si_force_per_length_to_lb(G);

	RBCParameters<double> rbc_params = create_rbc_params(
			shape,
			G_lb,
			unit_converter.si_force_per_length_to_lb(K),
			K_bend_rel * G_lb * shape->get_area());

	double a_eqv = std::sqrt(shape->get_area() / (4 * M_PI));

	shape->print_statistics();

	RBCParticle<double> rbc(shape_library.get_by_tag("RBC"), rbc_params);

	double vol_start = rbc.shape()->get_volume();
	double m = 0.1;				// Mass assigned to each node
	bool success = false;
	double vol_final = 0.59*vol_start;
	int N_it = 40000;

	plb::pcout << "Calibrating RBC shape" << std::endl;
	rbc.compute_forces();
	while(! success) {
		RBCParticle<double> rbc2 = rbc;
		success = true;

		int it = 0;
		bool converged = false;
		double Uref;
		while(!converged) {

			if(it <= N_it) {
				rbc2.params().vol_desired = vol_start + ((double)it / (double)N_it) * (vol_final - vol_start);
				rbc2.relax_nodes(m, 1.e-2);
			} else {
				rbc2.relax_nodes(m, 1.e-2);
			}

			// Check if the iteration has diverged (i.e. if area == NaN)
			// The c++ floating point specification defines that NaN != NaN.
			if(rbc2.area() != rbc2.area()) {
				m *= 2;
				success = false;
				plb::pcout << "Calibration failed, retrying with m = " << m << std::endl;
				break;
			}

			// Compute max velocity magnitude
			typedef RBCParticle<double>::vertex_const_iterator vertex_iterator;
			double max_vel_norm = 0;
			for(vertex_iterator iter = rbc2.begin(); iter != rbc2.end(); ++iter)
				max_vel_norm = std::max(max_vel_norm, std::sqrt(plb::norm(iter->vel)));

			if(it == N_it)
				Uref = max_vel_norm;
			else if(it > N_it && max_vel_norm/Uref < 0.1) {
				plb::pcout << "Calibration converged after " << it << " iteration" << std::endl;
				converged = true;
				break;
			}

			++it;

			if((it % 1000) == 0) {
				std::stringstream ss;
				ss << global::directories().getOutputDir() << "particleB_Kb" << K_bend_rel << "_" << it << ".vtu";
				rbc2.write_vtk(ss.str());
			}

			// Evaluate velocity field
			if((it % 1000) == 0)
				plb::pcout << "It: "<< it << ", total volume = " << rbc2.volume() << " (" << rbc2.params().vol_desired << ")"
						<< ", total area = " << rbc2.area() << " (" << rbc2.shape()->get_area() << ")"
						<< ", max velocity norm = " << max_vel_norm << std::endl;
		}

		if(success)
			rbc = rbc2;
	}

	/*shrink_rbc_volume(rbc, 0.59*rbc.volume(), 40000);*/

	rbc.set_minor_axis_orientation(Array<double,3>(0, 0, 1));
	rbc.set_center_of_mass(Array<double, 3>(0, 0, 0));

	std::stringstream ss;
	ss << global::directories().getOutputDir() << "particleB_Kb" << K_bend_rel << ".gmsh";
	rbc.write_gmsh(ss.str());

	/*pcout << "RBC Parameters:" << std::endl;
	pcout << "  Shear modulus: " << rbc_params.G() << std::endl;
	pcout << "  Area modulus " << rbc_params.K() << std::endl;
	pcout << "  Young's modulus: " << rbc_params.youngs_modulus() << std::endl;
	pcout << "  Poisson ratio " << rbc_params.poisson_ratio() << std::endl;
	pcout << "  Area dilation energy " << rbc_params.C << std::endl;
	pcout << "  Relative bending stiffness (kb/GA): " << rbc_params.k_bend / (rbc_params.G() * rbc.area()) << std::endl;
	pcout << "==============" << std::endl;
	pcout << "  Global area stiffness: " << rbc_params.k_area_global << std::endl;
	pcout << "  Local area stiffness: " << rbc_params.k_area_local << std::endl;
	pcout << "  Global volume stiffness: " << rbc_params.k_volume << std::endl;
	pcout << "  Bending stiffness: " << rbc_params.k_bend << std::endl;
	pcout << "  In-plane shear stiffness: " << rbc_params.k_in_plane << std::endl;
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
	pcout << "==============" << std::endl;*/

}
