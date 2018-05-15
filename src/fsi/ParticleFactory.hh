/*
 * ParticleFactory.hh
 *
 *  Created on: Jul 22, 2015
 *      Author: niber
 */

#ifndef PARTICLEFACTORY_HH_
#define PARTICLEFACTORY_HH_
#include "ParticleFactory.h"

namespace plb {

namespace fsi {

/*** ParticleCreator ***/
namespace detail {
template<class T, class ParticleType>
ParticleType * ParticleCreator<T, ParticleType>::create(const ParticleShape<T> * shape)
{
	return new ParticleType(shape);
}
}

/*** ParticleFactory ***/
template<class T>
ParticleFactory<T>::ParticleFactory()
: next_id(0), factories_()
{

}


template<class T>
plint ParticleFactory<T>::register_factory(detail::ParticleCreatorBase<T> * factory)
{
	factories_.push_back(factory);
	return factories_.size()-1;
}

template<class T>
void ParticleFactory<T>::clean_up()
{
	for(typename container_type::iterator it = factories_.begin(); it != factories_.end(); ++it)
		delete *it;
	factories_.clear();
}

/*template<class T>
	template<class BufferType>
ParticleBase3D<T> * ParticleFactory<T>::create(BufferType & buff, const ParticleShapeLibrary<T> & shape_library)
{
	// TODO implement an object register, so DeformableParticle does not need to know about all its subclasses
	plint type_id, shape_id;
	utils::unpack(buff, type_id);
	utils::unpack(buff, shape_id);

	ParticleBase3D<T> * ret = 0;
	const ParticleShape<T> * shape = shape_library.get_by_id(shape_id);

	typename container_type::iterator it = factories_.find(type_id);

	if(it == factories_.end()) {
		std::cerr << "ParticleBase3D<T>::create: unknown particle type " << type_id << std::endl;
	} else {
		it->second->create(shape);
		ret->unpack(buff);
	}

	return ret;
}*/

/*** Global methods ***/
template<class T>
ParticleFactory<T> & particleFactory()
{
	static ParticleFactory<T> instance;
	return instance;
}

template<class T, class ParticleType>
plint register_particle_type()
{
	return particleFactory<T>().register_factory(new detail::ParticleCreator<T, ParticleType>());
}


}

}



#endif /* PARTICLEFACTORY_HH_ */
