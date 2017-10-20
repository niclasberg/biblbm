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
			rbc2.relax_nodes(m, damping);

			// Check if the iteration has diverged (i.e. if area == NaN)
			// The c++ floating point specification defines that NaN != NaN.
			if(rbc2.area() != rbc2.area()) {
				m *= 2;
				success = false;
				plb::pcout << "Calibration failed, retrying with m = " << m << std::endl;
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

			++it;

			//if((it % 100) == 0)
			//	rbc2.print_energies(plb::pcout);

			if((it % 10000) == 0) {
				plb::pcout << "  It: "<< it << ", total volume = " << rbc2.volume() << " (" << rbc2.params().vol_desired << ")"
						<< ", total area = " << rbc2.area() << " (" << rbc2.shape()->get_area() << ")"
						<< ", max velocity norm = " << max_vel_norm << std::endl;
			}
		}

		if(success)
			rbc = rbc2;
	}
}

template<class T>
void create_rbc_minor_axis_transformation(plb::fsi::RBCParticle<T> & rbc, const plb::Array<T, 3> & pref_dir, plb::fsi::Transform<T> & transform)
{
	using namespace plb;
	using namespace fsi;


}


#endif /* CALIBRATE_RBC_SHAPE_H_ */
