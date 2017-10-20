/*
 * Particle.h
 *
 *  Created on: 22 feb 2014
 *      Author: niber
 */

#ifndef RIGIDPARTICLE_H_
#define RIGIDPARTICLE_H_
#include "core/globalDefs.h"
#include "core/array.h"
#include "core/geometry3D.h"
#include "ParticleShape.h"
#include "Quaternion.h"
#include "geometry.h"

namespace plb {

namespace fsi {

// Forward declarations
template<class T> class Boundary;

template<class T>
struct ParticleNodeInfo {
	Array<T, 3> normal;
	T area;
	Array<T, 3> pos;
	Array<T, 3> pos_rel;
	Array<T, 3> vel;
};

template<class T>
struct Transform {
public:
	Transform();
	Array<T, 3> apply(const Array<T, 3> &) const;
	void apply(const Array<T, 3> &, Array<T, 3> &) const;
	Transform & rotate(const Matrix<T, 3> &);
	Transform & scale(T);
	Transform combine_with(const Transform &);
	Transform & translate(const Array<T, 3> &);
//private:
	Matrix<T, 4> m;
};

template<class T>
class RigidParticle3D {
public:
	RigidParticle3D(const ParticleShape<T> * shape_)
	: shape(shape_), position(0, 0, 0), velocity(0, 0, 0),
	  ang_velocity(0, 0, 0), fsi_force(0, 0, 0), fsi_torque(0, 0, 0), density(1), id(-1),
	  coll_force(0, 0, 0), coll_force_last(0, 0, 0),
	  coll_torque(0, 0, 0), coll_torque_last(0, 0, 0),
	  main_processor(-1), ang_momentum(0, 0, 0), scale(1)
	{
		collision_forces.resize(shape->get_num_elements());
		for(int i = 0; i < collision_forces.size(); ++i)
			collision_forces[i].resetToZero();
		orientation.set_to_unity();
		orientation.to_rot_matrix(rot_matrix);
		update();
	}

	RigidParticle3D * clone() { return new RigidParticle3D(*this); }

	// Dynamic state setters
	void set_velocity(const Array<T, 3> & vel) { this->velocity = vel; }
	void set_angular_velocity(const Array<T, 3> & ang_vel) { this->ang_velocity = ang_vel; update(); }
	void set_position(const Array<T, 3> & pos) { this->position = pos; }
	void set_orientation_angle_axis(const T & angle, const Array<T, 3> & axis)
	{
		orientation.set_from_angle_axis(angle, axis); update();
	}

	void update_rotation_matrix() { orientation.to_rot_matrix(rot_matrix); }
	void update();

	// Kinetics
	void move(T damping = 0);
	void integrate_momentum_no_fsi(plint, T damping = 0);
	void integrate_momentum_1st_stage(plint, const Array<T, 3> &);
	void integrate_momentum_2nd_stage(plint);
	T compute_kinetic_energy() const;

	// Getters
	Array<T, 3> & get_position() { return position; }
	const Array<T, 3> & get_position() const { return position; }
	Array<T, 3> & get_velocity() { return velocity; }
	const Array<T, 3> & get_velocity() const { return velocity; }
	Array<T, 3> & get_angular_velocity() { return ang_velocity; }
	const Array<T, 3> & get_angular_velocity() const { return ang_velocity; }
	Quaternion<T> & get_orientation() { return orientation; }
	const Quaternion<T> & get_orientation() const { return orientation; }
	Array<T, 3> & get_force() { return fsi_force; }
	const Array<T, 3> & get_force() const { return fsi_force; }
	Array<T, 3> & get_torque() { return fsi_torque; }
	const Array<T, 3> & get_torque() const { return fsi_torque; }
	T & get_density() { return density; }
	const T & get_density() const { return density; }
	T & get_scale() { return scale; }
	const T & get_scale() const { return scale; }


	// Rigid body getters
	Box3D get_bounding_box() const;
	geo::Sphere<T> get_bounding_sphere() const { return geo::Sphere<T>(shape->get_radius(), position); }
	T get_radius() const { return shape->get_radius(); }
	pluint get_shape_id() const { PLB_PRECONDITION(shape); return shape->get_id(); }
	Matrix<T, 3> compute_inertia_world_frame() const;
	T get_mass() const { return std::pow(scale, 3) * shape->get_volume() * density; }
	Transform<T> get_fixed_to_world_transform() const;
	Transform<T> get_world_to_fixed_transform() const;

	// Node getters
	pluint get_vertex_count() const { return shape->get_num_elements(); }
	Array<T, 3> get_centroid_body_frame(pluint indx) const { return shape->get_centroid(indx); }
	void get_centroid_body_frame(pluint indx, Array<T, 3> & ret) const { ret = shape->get_centroid(indx); }
	Array<T, 3> get_centroid_world_frame(pluint idx) const { return position + rot_matrix * shape->get_centroid(idx); }
	void get_centroid_world_frame(pluint idx, Array<T, 3> & ret) const { ret = position + rot_matrix * shape->get_centroid(idx); }

	void get_node_info(pluint idx, ParticleNodeInfo<T> & ret) const;

	void get_node_velocity_world_frame(pluint indx, Array<T, 3> & ret) const
	{
		ret = velocity + crossProduct(ang_velocity, rot_matrix * shape->get_centroid(indx));
	}

	Array<T, 3> get_node_velocity_world_frame(pluint indx) const
	{
		Array<T, 3> ret; get_node_velocity_body_frame(indx, ret); return ret;
	}

	void get_normal(pluint idx, Array<T, 3> & ret) const { ret = rot_matrix * shape->get_node(idx).normal; }
	Array<T, 3> get_normal(pluint idx) const { return rot_matrix * shape->get_node(idx).normal; }
	T get_vertex_area(pluint idx) const { return shape->get_node(idx).area; }

	// Communication methods
	void set_id(plint id) { this->id = id; }
	plint get_id() const { return id; }
	void set_proc_id(plint id) { this->main_processor = id; }
	plint get_proc_id() const { return main_processor; }

	// Collision handling
	template<class InteractionPotential>
	bool compute_collision_forces(RigidParticle3D<T> &, const InteractionPotential &);
	template<class InteractionPotential, class Arithmetic>
	bool compute_collision_forces(RigidParticle3D<T> &, const InteractionPotential &, const Arithmetic &);
	template<class InteractionPotential>
	bool compute_wall_collision_forces(const Boundary<T> &, const InteractionPotential &);

	// IO
	void write_to_stream_as_vtk(std::ostream &) const;

	// temporary storage
	Array<T, 3> fsi_force_tmp;							// Temporary force (used in 2nd stage of fsi integration), world frame
	Array<T, 3> fsi_torque_tmp;							// Temporary torque, world frame
	Array<T, 3> coll_force, coll_force_last;		// Collision forces (world frame)
	Array<T, 3> coll_torque, coll_torque_last;		// Collision torques (world frame)
	Array<T, 3> ang_momentum;

	// For output
	std::vector<Array<T, 3> > collision_forces;
	std::vector<Array<T, 3> > hydrodynamic_forces;

private:
	void compute_angular_velocity();

	// Dynamic state
	Array<T, 3> position;			// Center of mass position (world frame)
	Array<T, 3> velocity;			// Center of mass velocity (world frame)
	Array<T, 3> ang_velocity;		// Angular velocity (world frame)
	fsi::Quaternion<T> orientation;	// Orientation quaternion
	Matrix<T, 3> rot_matrix;		// Rotation matrix (fixed to world frame)
	Array<T, 3> fsi_force;				// Force (world frame)
	Array<T, 3> fsi_torque;				// Torque (body frame)

	// Underlying shape
	const ParticleShape<T> * shape;	// Rigid body shape
	T density;						// Mass density
	T scale; 						// Rescaling of the underlying shape

	// Communication variables
	plint id;
	plint main_processor;
};

} /* namespace fsi */

} /* namespace plb */

#endif /* RIGIDPARTICLE_H_ */
