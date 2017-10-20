/*
 * Particle.h
 *
 *  Created on: 22 feb 2014
 *      Author: niber
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_
#include "core/globalDefs.h"
#include "core/array.h"

namespace plb {

template<class T>
class RigidParticle3D {
public:

	// Rigid body getters
	Box3D get_bounding_box() const;
	T get_radius() const;

	// Node getters
	pluint get_vertex_count() const;
	Array<T, 3> get_vertex_body_frame(pluint) const;
	void get_vertex_body_frame(pluint, Array<T, 3> &) const;
	Array<T, 3> get_vertex_world_frame(pluint) const;
	void get_vertex_world_frame(pluint, Array<T, 3> &) const;
	Array<T, 3> get_vertex_velocity_body_frame(pluint) const;
	void get_vertex_velocity_body_frame(pluint, Array<T, 3> &) const;
	Array<T, 3> get_vertex_velocity_world_frame(pluint) const;
	void get_vertex_velocity_world_frame(pluint, Array<T, 3>) const;
	T get_vertex_area(pluint);

private:
	// Dynamic state
	Array<T, 3> position; 		// Center of mass position (world frame)
	Array<T, 3> velocity;		// Center of mass velocity (world frame)
	Array<T, 3> ang_velocity;	// Angular velocity (world frame)
	void * orientation;			// Orientation quaternion (world frame)
	Array<T, 3> force;			// Force (world frame)
	Array<T, 3> torque;			// Torque (world frame)

	void * shape;				// Rigid body shape
};

}




#endif /* PARTICLE_H_ */
