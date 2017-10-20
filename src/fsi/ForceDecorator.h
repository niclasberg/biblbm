/*
 * ForceDecorator.h
 *
 *  Created on: Apr 12, 2016
 *      Author: niber
 */

#ifndef FORCEDECORATOR_H_
#define FORCEDECORATOR_H_
#include "ParticleBase.h"
#include "Boundary.h"

namespace plb {

namespace fsi {

template<class T>
struct ForceDecorator {
	virtual void apply_force(ParticleBase3D<T> * particle) = 0;
};

template<class T, class Interaction>
class WallInteraction : public ForceDecorator<T> {
public:
	WallInteraction(const Boundary<T> & boundary, const Interaction & interaction);
	virtual void apply_force(ParticleBase3D<T> * particle);
private:
	const Boundary<T> * boundary_;
	Interaction interaction_;
};

}

}




#endif /* FORCEDECORATOR_H_ */
