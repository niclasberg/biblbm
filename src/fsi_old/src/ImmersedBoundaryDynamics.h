/*
 * ImmersedBoundaryDynamics.h
 *
 *  Created on: 22 feb 2014
 *      Author: niber
 */

#ifndef IMMERSEDBOUNDARYDYNAMICS_H_
#define IMMERSEDBOUNDARYDYNAMICS_H_
#include "Particle.h"
#include <vector>

namespace plb {

template<class T, template<typename U> class Descriptor>
class ImmersedBoundaryDynamics3D {
public:
	void init();
	void compute_particle_motion();
	void compute_fsi_forces(const MultiBlockLattice3D<T, Descriptor> &, const MultiTensorField3D<T, 3> &);
	void compute_particle_interaction();
	void add_particle(Particle3D<T> *);
	void save_checkpoint();
	void load_checkpoint();

private:
	std::vector<Particle3D<T> *> particles;
};

}

#endif /* IMMERSEDBOUNDARYDYNAMICS_H_ */
