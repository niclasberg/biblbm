/*
 * ForceDecorator.hh
 *
 *  Created on: Apr 12, 2016
 *      Author: niber
 */

#ifndef FORCEDECORATOR_HH_
#define FORCEDECORATOR_HH_
#include "ForceDecorator.h"

namespace plb {

namespace fsi {

// Wall interaction force
template<class T, class Interaction>
WallInteraction<T, Interaction>::WallInteraction(const Boundary<T> & boundary, const Interaction & interaction)
: boundary_(&boundary), interaction_(interaction)
{

}

template<class T, class Interaction>
void WallInteraction<T, Interaction>::apply_force(ParticleBase3D<T> * particle)
{
	for(typename ParticleBase3D<T>::vertex_iterator it = particle->begin(); it != particle->end(); ++it) {
		T dist = boundary_->distance_to_boundary(it->pos);
		if(dist <= interaction_.get_cutoff_distance())
			it->force -= interaction_(dist) * (boundary_->get_normal(it->pos));
	}
}

}

}



#endif /* FORCEDECORATOR_HH_ */
