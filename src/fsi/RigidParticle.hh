#ifndef RIGIDPARTICLE_HH_
#define RIGIDPARTICLE_HH_
#include "RigidParticle.h"
#include "ParticleFactory.h"

namespace plb {

namespace fsi {

template<class T>
plint RigidParticle3D<T>::type_id = register_particle_type<T, RigidParticle3D<T> >();

template<class T>
RigidParticle3D<T>::RigidParticle3D(const ParticleShape<T> * shape)
: ParticleBase3D<T>(shape),
  velocity_(0, 0, 0),
  ang_velocity_(0, 0, 0),
  orientation_(),
  force_(0, 0, 0),
  torque_(0, 0, 0),
  density_(1),
  damping_(0.),
  relaxation_factor_(0.45),
  scale_(1.)
{
	orientation_.set_to_unity();
	update();
}

template<class T>
void RigidParticle3D<T>::update()
{
	// Generate new node positions
	update_rotation_matrix();
	for(plint i = 0; i < this->count_nodes(); ++i) {
		this->get_node(i).pos = this->center_of_mass() + rot_matrix_ * (scale_ * this->shape()->get_vertex(i));
	}
	ParticleBase3D<T>::update();
}

template<class T>
void RigidParticle3D<T>::update_rotation_matrix()
{
	orientation_.to_rot_matrix(this->rot_matrix_);
}

template<class T>
void RigidParticle3D<T>::reset_forces()
{
	ParticleBase3D<T>::reset_forces();
	this->force().resetToZero();
	this->torque().resetToZero();
}

template<class T>
void RigidParticle3D<T>::move_vertices(Boundary<T> * boundary)
{
	// Compute total torque and force
	this->force().resetToZero();
	this->torque().resetToZero();
	for(plint i = 0; i < this->count_nodes(); ++i) {
		const Vertex<T> & n = this->get_node(i);
		force_ += n.force;
		torque_ += crossProduct(n.pos - this->center_of_mass(), n.force);
	}

	const Array<T, 3> ang_velocity_body_frame = orientation_.apply_inv_rotation(ang_velocity_);
	const Matrix<T, 3> inertia_body_frame = (density_ * std::pow(scale_, 5)) * this->shape()->get_moment_of_inertia();
	//const Matrix<T, 3> inertia_world_frame = rot_matrix_.transpose() * inertia_body_frame * rot_matrix_;

	// Compute accelerations (world frame)
	const Array<T, 3> v_dot = force_ / get_mass() - damping_*velocity_;
	const Array<T, 3> omega_dot = orientation_.apply_rotation(
						inertia_body_frame.solve(
								orientation_.apply_inv_rotation(torque_) -
								crossProduct(ang_velocity_body_frame, inertia_body_frame * ang_velocity_body_frame)))
						- damping_*ang_velocity_;

	// Integrate linear and angular velocity
	velocity_ += v_dot;
	ang_velocity_ += omega_dot;

	// Integrate position
	this->center_of_mass_ += velocity_;

	// Quaternion integration. A unit norm-preserving exponential map is used.
	const Array<T, 3> ang_vel2 = ang_velocity_;
	const T om_norm = norm(ang_vel2);
	if(om_norm >= (T)10 * std::numeric_limits<T>::epsilon()) // avoid singularity when the norm of the angular velocity is small
		orientation_ = (Quaternion<T>(om_norm, ang_vel2 / om_norm)) * orientation_;

	// Update nodes
	this->update();
}

template<class T>
T RigidParticle3D<T>::get_mass() const
{
	return (scale_*scale_*scale_) * this->shape()->get_volume() * density_;
}

template<class T>
void RigidParticle3D<T>::compute_forces()
{
	// Compute weights
	std::vector<T> M(this->count_nodes(), T());
	for(typename ParticleShape<T>::triangle_const_iterator it = this->shape()->triangles_begin();
				it != this->shape()->triangles_end(); ++it) {
		T local_area = tri::triangle_area(this->get_node(it->i0).pos, this->get_node(it->i1).pos, this->get_node(it->i2).pos);
		M[it->i0] += local_area / 3.;
		M[it->i1] += local_area / 3.;
		M[it->i2] += local_area / 3.;
	}

	// Compute fsi forces
	for(plint i = 0; i < this->count_nodes(); ++i) {
		Vertex<T> & n = this->get_node(i);
		n.force += relaxation_factor_ * M[i] * (compute_node_velocity(i) - n.vel);
	}
}

template<class T>
Array<T, 3> RigidParticle3D<T>::compute_node_velocity(plint i) const
{
	return this->velocity() +
			crossProduct(this->angular_velocity(),
					this->get_node(i).pos - this->center_of_mass());
}

template<class T>
void RigidParticle3D<T>::pack(std::vector<char> & buf) const
{
	ParticleBase3D<T>::pack(buf);

	utils::pack(buf, this->center_of_mass_);
	utils::pack(buf, velocity_);
	utils::pack(buf, ang_velocity_);
	utils::pack(buf, orientation_);
	utils::pack(buf, force_);
	utils::pack(buf, torque_);
	utils::pack(buf, density_);
	utils::pack(buf, damping_);
	utils::pack(buf, relaxation_factor_);
	utils::pack(buf, scale_);
}

template<class T>
void RigidParticle3D<T>::unpack(char *& buf)
{
	ParticleBase3D<T>::unpack(buf);

	utils::unpack(buf, this->center_of_mass_);
	utils::unpack(buf, velocity_);
	utils::unpack(buf, ang_velocity_);
	utils::unpack(buf, orientation_);
	utils::unpack(buf, force_);
	utils::unpack(buf, torque_);
	utils::unpack(buf, density_);
	utils::unpack(buf, damping_);
	utils::unpack(buf, relaxation_factor_);
	utils::unpack(buf, scale_);

	this->update();
}

template<class T>
void RigidParticle3D<T>::unpack(std::istream & buf)
{
	ParticleBase3D<T>::unpack(buf);

	utils::unpack(buf, this->center_of_mass_);
	utils::unpack(buf, velocity_);
	utils::unpack(buf, ang_velocity_);
	utils::unpack(buf, orientation_);
	utils::unpack(buf, force_);
	utils::unpack(buf, torque_);
	utils::unpack(buf, density_);
	utils::unpack(buf, damping_);
	utils::unpack(buf, relaxation_factor_);
	utils::unpack(buf, scale_);

	this->update();
}

} /* namespace fsi */

} /* namespace plb */

#endif /* RIGIDPARTICLE_HH_ */
