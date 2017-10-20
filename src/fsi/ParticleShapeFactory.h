/*
 * ParticleShapeFactory.h
 *
 *  Created on: Mar 4, 2014
 *      Author: niber
 */

#ifndef PARTICLESHAPEFACTORY_H_
#define PARTICLESHAPEFACTORY_H_
#include "ParticleShape.h"
#include "core/globalDefs.h"
#include <string>
#include <map>

namespace plb {

namespace fsi {

template<class T>
class ParticleShapeLibrary {
public:
	ParticleShapeLibrary() : next_id(0) { }
	~ParticleShapeLibrary();

	void clear();
	void store_mesh(Array<T, 3> *, unsigned int *, unsigned int, unsigned int, const std::string & tag);
	void read_and_store_mesh(const std::string & file_name, const std::string & tag);
	ParticleShape<T> * get_by_tag(std::string tag);
	const ParticleShape<T> * get_by_tag(std::string tag) const;
	const ParticleShape<T> * get_by_id(plint id) const;
	ParticleShape<T> * get_by_id(plint id);
	plint get_id_by_tag(std::string tag) const;
	plint get_shape_count() const;
	T get_max_particle_radius() const;

private:
	plint next_id;
	std::vector<ParticleShape<T> *> id_shape_map;
	std::map<std::string, plint> name_id_map;
};

}

}




#endif /* PARTICLESHAPEFACTORY_H_ */
