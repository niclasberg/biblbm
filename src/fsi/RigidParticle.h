/*
 * Particle.h
 *
 *  Created on: 22 feb 2014
 *      Author: niber
 */

#ifndef RIGIDPARTICLE_H_
#define RIGIDPARTICLE_H_
#include "ParticleBase.h"

namespace plb {

namespace fsi {

// Forward declarations
template<class T> class Boundary;

template<class T>
class RigidParticle3D : public ParticleBase3D<T> {
public:
	typedef typename ParticleBase3D<T>::vertex_const_iterator vertex_const_iterator;
	typedef typename ParticleBase3D<T>::vertex_iterator vertex_iterator;

	RigidParticle3D(const ParticleShape<T> * shape);
	virtual ~RigidParticle3D() { }
	virtual RigidParticle3D * clone() const { return new RigidParticle3D(*this); }

	virtual void update();
	virtual void move_vertices(Boundary<T> * boundary = 0);
	virtual void compute_forces();
	virtual void reset_forces();

	// Getters
	Array<T, 3> & velocity() { return velocity_; }
	const Array<T, 3> & velocity() const { return velocity_; }
	Array<T, 3> & angular_velocity() { return ang_velocity_; }
	const Array<T, 3> & angular_velocity() const { return ang_velocity_; }
	Quaternion<T> & orientation() { return orientation_; }
	const Quaternion<T> & orientation() const { return orientation_; }
	Array<T, 3> & force() { return force_; }
	const Array<T, 3> & force() const { return force_; }
	Array<T, 3> & torque() { return torque_; }
	const Array<T, 3> & torque() const { return torque_; }
	T & density() { return density_; }
	const T & density() const { return density_; }
	T & scale() { return scale_; }
	const T & scale() const { return scale_; }
	T & damping() { return damping_; }
	const T & damping() const { return damping_; }

	// Rigid body getters
	Transform<T> get_fixed_to_world_transform() const;
	Transform<T> get_world_to_fixed_transform() const;

	// IO
	virtual plint get_type_id() const { return type_id; }
	virtual void pack(std::vector<char> &) const;
	virtual void unpack(char *&);
	virtual void unpack(std::istream &);
	virtual void write_lightweight(std::ostream &) const { }

private:
	Array<T, 3> compute_node_velocity(plint i) const;
	Array<T, 3> compute_node_position(plint i) const;
	void update_rotation_matrix();
	T get_mass() const;

private:
	static plint type_id;

	// Dynamic state
	Array<T, 3> velocity_;				// Center of mass velocity (world frame)
	Array<T, 3> ang_velocity_;			// Angular velocity (world frame)
	fsi::Quaternion<T> orientation_;	// Orientation quaternion
	Matrix<T, 3> rot_matrix_;			// Rotation matrix (fixed to world frame)
	Array<T, 3> force_;					// Force (world frame)
	Array<T, 3> torque_;				// Torque (body frame)
	T density_;							// Mass density
	T damping_;
	T relaxation_factor_;
	T scale_;
};

} /* namespace fsi */

} /* namespace plb */

#endif /* RIGIDPARTICLE_H_ */
