/*
 * calibrate_rbc_shape.h
 *
 *  Created on: Jun 30, 2015
 *      Author: niber
 */

#ifndef CALIBRATE_RBC_SHAPE_H_
#define CALIBRATE_RBC_SHAPE_H_
#include "fsi/RBCParticle.h"

template<class T>
void shrink_rbc_volume(plb::fsi::RBCParticle<T> & rbc, T vol_final, int N_it, T m = 1, T damping = 1.e-2)
{
	typedef typename plb::fsi::RBCParticle<T>::vertex_const_iterator vertex_iterator;

	T vol_start = rbc.shape()->get_volume();
	bool success = false;

	plb::pcout << "Calibrating RBC shape" << std::endl;
	rbc.compute_forces();
	while(! success) {
		plb::fsi::RBCParticle<T> rbc2 = rbc;
		success = true;

		int it = 0;
		bool converged = false;
		T Uref;
		while(!converged) {
			if(it <= N_it)
				rbc2.params().vol_desired = vol_start + ((T)it / (T)N_it) * (vol_final - vol_start);
			rbc2.relax_nodes(damping);

			// Check if the iteration has diverged (i.e. if area == NaN)
			// The c++ floating point specification defines that NaN != NaN.
			if(rbc2.area() != rbc2.area()) {
				damping *= 2.;
				success = false;
				plb::pcout << "Calibration failed, retrying with damping = " << damping << std::endl;
				break;
			}

			// Compute max velocity magnitude
			T max_vel_norm = 0;
			for(vertex_iterator iter = rbc2.begin(); iter != rbc2.end(); ++iter)
				max_vel_norm = std::max(max_vel_norm, std::sqrt(plb::norm(iter->vel)));

			/*if(it == N_it)
				Uref = max_vel_norm;
			else if(it > N_it && max_vel_norm/Uref < 1.e-2) {
				plb::pcout << "Calibration converged after " << it << " iteration" << std::endl;
				converged = true;
				break;
			}*/
			if(it == 2*N_it)
				converged = true;

			/*if((it % 100) == 0)
				rbc2.print_energies(plb::pcout);*/

			if((it % 10000) == 0) {
				plb::pcout << "  It: "<< it << ", total volume = " << rbc2.volume() << " (" << rbc2.params().vol_desired << ")"
						<< ", total area = " << rbc2.area() << " (" << rbc2.shape()->get_area() << ")"
						<< ", max velocity norm = " << max_vel_norm << std::endl;
			}
			++it;
		}

		if(success)
			rbc = rbc2;
	}
}

template<class T>
plb::fsi::RBCParameters<T> create_rbc_params(
	const plb::fsi::ParticleShape<T> * shape,
	T shear_modulus,
	T K_area,
	T K_bend
)
{
	using namespace plb;
	using namespace fsi;

	RBCParameters<T> rbc_params;

	// Get average link length
	T lavg = 0;
	for(typename ParticleShape<T>::link_const_iterator it = shape->links_begin(); it != shape->links_end(); ++it)
		lavg += it->length;
	lavg /= shape->count_links();

	// In plane properties
	rbc_params.shear_modulus = shear_modulus;
	rbc_params.L0 = lavg;
	rbc_params.vol_desired = shape->get_volume();

	// Area constraints
	rbc_params.k_area_global = 0;
	rbc_params.k_area_local = 0;
	T delta_k_area = K_area - rbc_params.K();

	rbc_params.k_area_global = delta_k_area * 0.8;
	rbc_params.k_area_local = delta_k_area * 0.2;
	rbc_params.k_volume = 1000. * shear_modulus / shape->get_radius();

	// Bending stiffness
	rbc_params.k_bend = K_bend;
	rbc_params.theta0 = std::acos((std::sqrt(3.)*(shape->count_vertices() - 2.) - 5*M_PI) / (std::sqrt(3.)*(shape->count_vertices() - 2.) - 3.*M_PI) );

	return rbc_params;
}

template<class T>
plb::fsi::RBCParameters<T> create_platelet_params(
	const plb::fsi::ParticleShape<T> * shape,
	T shear_modulus,
	T K_area,
	T K_bend
)
{
	using namespace plb;
	using namespace fsi;

	RBCParameters<T> rbc_params;

	// Get average link length
	T lavg = 0;
	for(typename ParticleShape<T>::link_const_iterator it = shape->links_begin(); it != shape->links_end(); ++it)
		lavg += it->length;
	lavg /= shape->count_links();

	// In plane properties
	rbc_params.shear_modulus = shear_modulus;
	rbc_params.L0 = lavg;
	rbc_params.vol_desired = shape->get_volume();

	// Area constraints
	//rbc_params.k_area_global = 0;
	//rbc_params.k_area_local = 0;
	//T delta_k_area = K_area - rbc_params.K();

	rbc_params.k_area_global = 0.1;
	rbc_params.k_area_local = 0.1;
	rbc_params.k_volume = 0.5;

	// Bending stiffness
	rbc_params.k_bend = K_bend;
	rbc_params.theta0 = std::acos((std::sqrt(3.)*(shape->count_vertices() - 2.) - 5*M_PI) / (std::sqrt(3.)*(shape->count_vertices() - 2.) - 3.*M_PI) );

	return rbc_params;
}

template<class T>
class RBCShapeInitializer {
public:
	RBCShapeInitializer(const std::string & mesh_file) : rbc(0), shape_library() 
	{ 
		shape_library.read_and_store_mesh(mesh_file, "RBC");
		rbc = new plb::fsi::RBCParticle<T>(shape_library.get_by_tag("RBC"));
	}

	RBCShapeInitializer(plb::fsi::RBCParticle<T> * rbc_) : rbc(rbc_), shape_library() 
	{ 

	}

	~RBCShapeInitializer() {
		delete rbc;
	}

	const plb::fsi::RBCParticle<T> & get_rbc() const { return *rbc; }
	plb::fsi::RBCParticle<T> & get_rbc() { return *rbc; }

	plb::fsi::RBCParameters<T> & get_params() { return rbc->params(); }
	const plb::fsi::RBCParameters<T> & get_params() const { return rbc->params(); }
	const plb::fsi::ParticleShape<T> * get_shape() const { return rbc->shape(); }

	void compute_rbc_parameters(T shear_modulus, T K_area, T K_bend) 
	{
		rbc->params() = create_rbc_params(rbc->shape(), shear_modulus, K_area, K_bend);
	}

	void shrink_volume(T vol_frac, plb::plint Nit) 
	{
		shrink_rbc_volume(*rbc, vol_frac * get_shape()->get_volume(), Nit);
		rbc->set_minor_axis_orientation(plb::Array<T, 3>(0, 0, 1));
		rbc->params().vol_desired = rbc->volume();
	}

	void store_shape(plb::fsi::ParticleShapeLibrary<T> & shape_lib, const char * tag) const
	{
		rbc->store_shape(shape_lib, tag);
	}
	
private:
	// Disable copy constructor and assignment operator
	RBCShapeInitializer & operator=(const RBCShapeInitializer &);
	RBCShapeInitializer(const RBCShapeInitializer &);

	plb::fsi::RBCParticle<T> * rbc;
	plb::fsi::ParticleShapeLibrary<T> shape_library; //Contains only the current rbc shape
};

/*template<class T>
void create_rbc_minor_axis_transformation(plb::fsi::RBCParticle<T> & rbc, const plb::Array<T, 3> & pref_dir, plb::fsi::Transform<T> & transform)
{
	using namespace plb;
	using namespace fsi;


}*/


#endif /* CALIBRATE_RBC_SHAPE_H_ */
