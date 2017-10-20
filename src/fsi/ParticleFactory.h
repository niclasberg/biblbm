/*
 * ParticleFactory.h
 *
 *  Created on: Jun 12, 2015
 *      Author: niber
 */

#ifndef PARTICLEFACTORY_H_
#define PARTICLEFACTORY_H_

#include "ParticleBase.h"
#include "ParticleShapeFactory.h"
#include "utils.h"

namespace plb {

namespace fsi {

template<class T> class ParticleBase3D;
template<class T> class ParticleShapeLibrary;

namespace detail {
template<class T>
struct ParticleCreatorBase {
	virtual ParticleBase3D<T> * create(const ParticleShape<T> * shape) = 0;
};

template<class T, class ParticleType>
struct ParticleCreator : public ParticleCreatorBase<T> {
	virtual ParticleType * create(const ParticleShape<T> * shape);
};

}

template<class T>
class ParticleFactory {
public:
	typedef std::vector<detail::ParticleCreatorBase<T> * > container_type;

	ParticleFactory();

	// Particle registration and creation
	template<class BufferType>
	ParticleBase3D<T> * create(BufferType & buff, const ParticleShapeLibrary<T> & shape_library)
	{
		plint type_id, shape_id;
		utils::unpack(buff, type_id);
		utils::unpack(buff, shape_id);

		ParticleBase3D<T> * ret = 0;
		const ParticleShape<T> * shape = shape_library.get_by_id(shape_id);

		if(type_id >= factories_.size() || type_id < 0) {
			std::cerr << "ParticleBase3D<T>::create: unknown particle type " << type_id << std::endl;
		} else {
			ret = factories_[type_id]->create(shape);
			ret->unpack(buff);
		}

		return ret;
	}

	plint register_factory(detail::ParticleCreatorBase<T> *);

public:
	container_type factories_;
	plint next_id;
};

template<class T>
ParticleFactory<T> & particleFactory();

template<class T, class ParticleType>
plint register_particle_type();

}

}



#endif /* PARTICLEFACTORY_H_ */
