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
void RigidParticle3D<T>::move_vertices()
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
	// Compute fsi forces
	for(plint i = 0; i < this->count_nodes(); ++i) {
		Vertex<T> & n = this->get_node(i);
		n.force += relaxation_factor_ * (n.vel - compute_node_velocity(i));
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

/*

template<class T>
void RigidParticle3D<T>::move(T damping)
{
	const Array<T, 3> ang_velocity_body_frame = orientation.apply_inv_rotation(ang_velocity);
	const Matrix<T, 3> & inertia = (density * std::pow(scale, 5)) * shape->get_moment_of_inertia();

	// Compute accelerations (world frame)
	const Array<T, 3> v_dot = (fsi_force + coll_force) / get_mass() - damping*velocity;
	const Array<T, 3> omega_dot = orientation.apply_rotation(
						inertia.solve(
								orientation.apply_inv_rotation(fsi_torque + coll_torque - damping*ang_momentum) -
								crossProduct(ang_velocity_body_frame, inertia * ang_velocity_body_frame)));

	// Integrate velocity and position
	position += velocity + 0.5 * v_dot;

	// Quaternion integration. A unit norm-preserving exponential map is used.
	const Array<T, 3> ang_vel2 = ang_velocity + 0.5 * omega_dot;
	const T om_norm = norm(ang_vel2);
	if(om_norm >= 1.e-12) // avoid singularity when the norm of the angular velocity is small
		orientation = (Quaternion<T>(om_norm, ang_vel2 / om_norm)) * orientation;

	// Quaternion integration using a corrected taylor series expansion
	//const Quaternion<T> omega = Quaternion<T>(0, ang_velocity[0], ang_velocity[1], ang_velocity[2]);
	//const Quaternion<T> omega_dot2 = Quaternion<T>(0, omega_dot[0], omega_dot[1], omega_dot[2]);
	//const Quaternion<T> q_dot = 0.5 * omega;
	//const Quaternion<T> q_dot_dot = 0.5 * (omega_dot2 * orientation + omega * q_dot);
	//const T correction = 1 - 0.5*q_dot.norm_sqr() - std::sqrt(1 - q_dot.norm_sqr() - q_dot.dot(q_dot_dot) - 0.25*(q_dot_dot.norm_sqr() - util::sqr(q_dot.norm_sqr())));
	//orientation = orientation + q_dot + 0.5*q_dot_dot - correction*orientation;

	orientation.to_rot_matrix(rot_matrix);
}

template<class T>
void RigidParticle3D<T>::integrate_momentum_no_fsi(plint it, T damping)
{
	// Integrate velocity
	velocity += (0.5 * (coll_force + coll_force_last)) / get_mass() - damping*velocity;

	// Integrate angular velocity
	ang_momentum += 0.5 * (coll_torque + coll_torque_last) - damping * ang_momentum;
	compute_angular_velocity();
}

template<class T>
void RigidParticle3D<T>::compute_angular_velocity()
{
	ang_velocity = compute_inertia_world_frame().solve(ang_momentum);
}

template<class T>
void RigidParticle3D<T>::integrate_momentum_1st_stage(plint it, const Array<T, 3> & external_force)
{
	velocity += 0.5*(fsi_force + coll_force + coll_force_last) / get_mass()
									- external_force / density;
	ang_momentum += 0.5 * (fsi_torque + coll_torque + coll_torque_last);
	compute_angular_velocity();
}

template<class T>
void RigidParticle3D<T>::integrate_momentum_2nd_stage(plint it)
{
	velocity += 0.5 * fsi_force_tmp / get_mass();
	ang_momentum += 0.5 * fsi_torque_tmp;
	compute_angular_velocity();
}

template<class T>
	template<class InteractionPotential>
bool RigidParticle3D<T>::compute_collision_forces(RigidParticle3D<T> & p2, const InteractionPotential & potential)
{
	return compute_collision_forces(p2, potential, NormalArithmetic<T>());
}

template<class T>
	template<class InteractionPotential, class Arithmetic>
bool RigidParticle3D<T>::compute_collision_forces(RigidParticle3D<T> & p2, const InteractionPotential & potential, const Arithmetic & a)
{
	// Center to center distance (world frame)
	Array<T, 3> dr = a.vec_diff(p2.position, this->position);
	T dr_norm_sqr = normSqr(dr);

	// Do a bounding sphere overlap test to reduce computational overhead
	if(dr_norm_sqr > util::sqr(get_scale()*get_radius() + p2.get_scale()*p2.get_radius() + potential.get_cutoff_distance()))
		return false;

	bool ret = false;

	// Unit vector between the center of masses
	Array<T, 3> dr_unit = dr / std::sqrt(dr_norm_sqr);

	// Compute the transformation matrix from this particle's fixed coordinate system to that of p2
	Transform<T> transform;
	transform.scale(scale)
			.rotate(rot_matrix)
			.translate(-dr)
			.rotate(p2.rot_matrix.transpose())
			.scale((T)1.0 / p2.scale);

	// Compute the rotation matrix p1 -> p2 (no translation)
	Matrix<T, 3> rotation = p2.rot_matrix.transpose() * this->rot_matrix;

	// The center to center unit vector in p1 and p2's coordinate system
	Array<T, 3> dr_unit_p1 = orientation.apply_inv_rotation(dr_unit);
	Array<T, 3> dr_unit_p2 = p2.orientation.apply_inv_rotation(dr_unit);

	Array<T, 3> pos, d;
	plint match;
	T scaled_cutoff_distance = potential.get_cutoff_distance() / p2.scale;

	for(pluint i = 0; i < shape->get_num_elements(); ++i) {
		// Omit elements facing in the direction opposite of dr
		if(dot(dr_unit_p1, shape->get_node(i).normal) < 0)
			continue;

		// Evaluate the position in p2's fixed frame
		transform.apply(shape->get_centroid(i), pos);

		// Evaluate the distance between the current node and the closest node on the other surface.
		const plint match = p2.shape->get_closest_element_facing_in_direction(pos, -dr_unit_p2, scaled_cutoff_distance);
		if(match != -1) {
			T dist = p2.scale * dot(p2.shape->get_centroid(match) - pos, dr_unit_p2);

			if(std::abs(dist) < potential.get_cutoff_distance()) {
				const T projected_area = shape->get_node(i).area *
									dot(dr_unit, get_normal(i)) *
									(-dot(dr_unit, p2.get_normal(match)));
				Array<T, 3> df = -potential(dist) * projected_area * dr_unit;

				collision_forces[i] += df / shape->get_node(i).area;
				p2.collision_forces[match] += df / p2.shape->get_node(match).area;

				this->coll_force += df;
				this->coll_torque += crossProduct(scale * rot_matrix * shape->get_centroid(i), df);
				p2.coll_force -= df;
				p2.coll_torque -= crossProduct(p2.scale * p2.rot_matrix * p2.shape->get_centroid(match), df);

				ret = true;
			}
		}
	}

	return ret;
}

template<class T>
	template<class InteractionPotential>
bool RigidParticle3D<T>::compute_wall_collision_forces(
		const Boundary<T> & boundary,
		const InteractionPotential & potential)
{
	// Broad phase test:
	if( ! boundary.distance_to_boundary_less_than(position, (scale * shape->get_radius()) + potential.get_cutoff_distance()))
		return false;

	// Get normal direction
	const Array<T, 3> h = boundary.get_normal(position);

	bool ret = false;
	Transform<T> transform = get_fixed_to_world_transform();
	Array<T, 3> pos, normal;

	for(plint i = 0; i < shape->get_num_elements(); ++i) {
		const ParticleShapeNode<T> & node = shape->get_node(i);
		normal = rot_matrix * node.normal;
		if(dot(normal, h) > 0) {
			transform.apply(node.centroid, pos);

			T dist = boundary.distance_to_boundary(pos, h);

			if(dist <= potential.get_repulsion_distance()) {
				const Array<T, 3> df = -h * dot(h, normal) * node.area * potential(dist);
				collision_forces[i] += df / node.area;
				coll_force += df;
				coll_torque += crossProduct(pos - position, df);
				ret = true;
			}
		}
	}

	return ret;
}

template<class T>
Transform<T> RigidParticle3D<T>::get_fixed_to_world_transform() const
{
	Transform<T> ret;
	return ret.scale(scale).rotate(rot_matrix).translate(position);
}

template<class T>
Transform<T> RigidParticle3D<T>::get_world_to_fixed_transform() const
{
	Transform<T> ret;
	return ret.translate(-position).rotate(rot_matrix.transpose()).scale(1 / scale);
}

template<class T>
void RigidParticle3D<T>::reset_collision_forces()
{
	for(plint i = 0; i < collision_forces.size(); ++i)
		collision_forces[i].resetToZero();
	coll_force.resetToZero();
	coll_torque.resetToZero();
}

template<class T>
T RigidParticle3D<T>::compute_minimal_distance(const RigidParticle3D<T> & p2) const
{
	return compute_minimal_distance(p2, NormalArithmetic<T>());
}

template<class T>
	template<class Arithmetic>
T RigidParticle3D<T>::compute_minimal_distance(const RigidParticle3D<T> & p2, const Arithmetic & a) const
{
	// Compute the transformation matrix from this particle's fixed coordinate system to that of p2
	Transform<T> transform;
	transform.scale(scale)
			.rotate(rot_matrix)
			.translate(-a.vec_diff(p2.position, this->position))
			.rotate(p2.rot_matrix.transpose())
			.scale((T)1.0 / p2.scale);
	Array<T, 3> diff, diff_min, r;
	T diff_mag_sqr;
	T diff_min_mag_sqr = std::numeric_limits<T>::max();

	// This method is really slow (O(N^2))
	for(plint i = 0; i < shape->get_num_elements(); ++i) {
		transform.apply(shape->get_centroid(i), r);
		for(plint j = 0; j < p2.shape->get_num_elements(); ++j) {
			diff = r - p2.shape->get_centroid(j);
			diff_mag_sqr = normSqr(diff);
			if(diff_mag_sqr < diff_min_mag_sqr) {
				diff_min_mag_sqr = diff_mag_sqr;
				diff_min = diff;
			}
		}
	}
	return std::sqrt(diff_min_mag_sqr);
}

*/
} /* namespace fsi */

} /* namespace plb */

#endif /* RIGIDPARTICLE_HH_ */
