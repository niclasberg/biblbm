/*
 * ParticleParticleInteraction.h
 *
 *  Created on: Apr 19, 2016
 *      Author: niber
 */

#ifndef PARTICLEPARTICLEINTERACTION_H_
#define PARTICLEPARTICLEINTERACTION_H_


namespace plb {

namespace fsi {

template<class T>
struct PPInteractionForce {
	T operator()(T dist) const { return eval(dist); }
	virtual T eval(T dist) const = 0;
	virtual T get_cutoff_distance() const = 0;
};

template<class T, class Interaction>
class PPPotentialForce : public PPInteractionForce<T> {
public:
	PPPotentialForce(const Interaction & potential) : potential_(potential) { }
	virtual T eval(T dist) const
	{
		return potential_(dist);
	}

	virtual T get_cutoff_distance() const
	{
		return potential_.get_cutoff_distance();
	}

private:
	Interaction potential_;
};

} /* namespace fsi */

} /* namespace plb */


#endif /* PARTICLEPARTICLEINTERACTION_H_ */
