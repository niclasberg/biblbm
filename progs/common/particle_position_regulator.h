#ifndef PARTICLE_POSITION_REGULATOR_H_
#define PARTICLE_POSITION_REGULATOR_H_

template<class T>
struct ParticlePositionRegulator : public plb::fsi::ForceDecorator<T> 
{
	ParticlePositionRegulator(const plb::Array<T, 3> & pos, T p_) : p(p_), position(pos) { }

	virtual void apply_force(plb::fsi::ParticleBase3D<T> * particle)
	{
		plb::Array<T, 3> force = p * (position - particle->center_of_mass());
		for(typename plb::fsi::ParticleBase3D<T>::vertex_iterator it = particle->begin(); it != particle->end(); ++it) 
			it->force += force / (T) particle->count_nodes();
	}

	T p;
	plb::Array<T, 3> position;
};

#endif /* PARTICLE_POSITION_REGULATOR_H_ */
