/*
 * rbc_parameters.h
 *
 *  Created on: Feb 2, 2016
 *      Author: niber
 */

#ifndef RBC_PARAMETERS_H_
#define RBC_PARAMETERS_H_
#include "fsi/ParticleShape.h"
#include "fsi/RBCParticle.h"

template<class T>
plb::fsi::RBCParameters<T> create_rbc_params(
	plb::fsi::ParticleShape<T> * shape,
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
	rbc_params.theta0 = std::acos((std::sqrt(3)*(shape->count_vertices() - 2) - 5*M_PI) / (std::sqrt(3)*(shape->count_vertices() - 2) - 3*M_PI) );

	return rbc_params;
}



#endif /* RBC_PARAMETERS_H_ */
