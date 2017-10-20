/*
 * ParticlePositionInitializer.hh
 *
 *  Created on: Apr 2, 2014
 *      Author: niber
 */

#ifndef PARTICLEPOSITIONINITIALIZER_HH_
#define PARTICLEPOSITIONINITIALIZER_HH_
#include "ParticlePositionInitializer.h"
#include <cstdlib>

namespace plb {

namespace fsi {

template<class T>
bool ParticlePositionInitializer<T>::generate_points(pluint N, T radius, char tag)
{
	if(2*radius > grid.get_min_cell_size()) {
		pcerr << "ParticlePositionGenerator:generate_points: particle radius is larger than the grid spacing!" << std::endl;
		return false;
	}

	ParticlePosition<T> pp;
	pp.tag = tag;
	pp.radius = radius;
	for(pluint it = 0; it < N; ++it) {
		bool success = false;
		for(plint i = 0; i < 5e6; ++i) {
			success = generate_position(pp);
			if(success)
				break;
		}

		if( ! success) {
			return false;
		}

		positions.push_back(pp);
		grid.insert(pp, pp.pos);
	}
	return true;
}

template<class T>
bool ParticlePositionInitializer<T>::generate_position(ParticlePosition<T> & pp)
{
	// Randomize position
	pp.pos[0] = (T) bounding_box.x0 + pp.radius + ((T) rand() / (T) RAND_MAX) * (bounding_box.x1 - bounding_box.x0 - 2*pp.radius);
	pp.pos[1] = (T) bounding_box.y0 + pp.radius + ((T) rand() / (T) RAND_MAX) * (bounding_box.y1 - bounding_box.y0 - 2*pp.radius);
	pp.pos[2] = (T) bounding_box.z0 + pp.radius + ((T) rand() / (T) RAND_MAX) * (bounding_box.z1 - bounding_box.z0 - 2*pp.radius);

	// Check if within boundary region
	if(boundary) {
		if( ! boundary->contains(pp.pos, pp.radius))
			return false;
	}

	Dot3D ci = grid.get_index(pp.pos);
	for(plint dx = -1; dx <= 1; ++dx) {
		for(plint dy = -1; dy <= 1; ++dy) {
			for(plint dz = -1; dz <= 1; ++dz) {
				Cell<T, ParticlePosition<T> > & cell = grid.get_cell(ci.x+dx, ci.y+dy, ci.z+dz);

				for(plint i = 0; i < cell.nodes.size(); ++i) {
					if(normSqr(cell.nodes[i].pos - pp.pos) < util::sqr(cell.nodes[i].radius + pp.radius)) {
						return false;
					}
				}
			}
		}
	}
	return true;
}

template<class T>
void ParticlePositionInitializer<T>::get_points(char tag, std::vector<Array<T, 3> > & points)
{
	for(typename std::vector<ParticlePosition<T> >::iterator it = positions.begin(); it != positions.end(); ++it) {
		if(it->tag == tag)
			points.push_back(it->pos);
	}
}

}

}



#endif /* PARTICLEPOSITIONINITIALIZER_HH_ */
