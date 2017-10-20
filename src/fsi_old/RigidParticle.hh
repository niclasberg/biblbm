#ifndef RIGIDPARTICLE_HH_
#define RIGIDPARTICLE_HH_
#include "RigidParticle.h"
#include "Boundary.h"

namespace plb {

namespace fsi {

/* Transform object methods */
template<class T>
Transform<T>::Transform()
{
	// Initialize to identity
	m(0, 0) = 1; m(0, 1) = 0; m(0, 2) = 0; m(0, 3) = 0;
	m(1, 0) = 0; m(1, 1) = 1; m(1, 2) = 0; m(1, 3) = 0;
	m(2, 0) = 0; m(2, 1) = 0; m(2, 2) = 1; m(2, 3) = 0;
	m(3, 0) = 0; m(3, 1) = 0; m(3, 2) = 0; m(3, 3) = 1;
}

template<class T>
Array<T, 3> Transform<T>::apply(const Array<T, 3> & v) const
{
	Array<T, 3> ret;
	apply(v, ret);
	return ret;
}

template<class T>
void Transform<T>::apply(const Array<T, 3> & v, Array<T, 3> & ret) const
{
	ret[0] = v[0]*m(0, 0) + v[1]*m(0, 1) + v[2]*m(0, 2) + m(0, 3);
	ret[1] = v[0]*m(1, 0) + v[1]*m(1, 1) + v[2]*m(1, 2) + m(1, 3);
	ret[2] = v[0]*m(2, 0) + v[1]*m(2, 1) + v[2]*m(2, 2) + m(2, 3);
}

template<class T>
Transform<T> & Transform<T>::rotate(const Matrix<T, 3> & r)
{
	// Equivalent to [r00 r01 r02 0 ; r10 r11 r12 0 ; r20 r21 r22 0 ; 0 0 0 1] * m
	T tmp[3];
	for(int j = 0; j < 4; ++j) {
		for(int i = 0; i < 3; ++i)
			tmp[i] = r(i, 0)*m(0, j) + r(i, 1)*m(1, j) + r(i, 2)*m(2, j);

		for(int i = 0; i < 3; ++i)
			m(i, j) = tmp[i];
	}

	return *this;
}

template<class T>
Transform<T> & Transform<T>::scale(T s)
{
	// Equivalent to [s 0 0 0 ; 0 s 0 0 ; 0 0 s 0 ; 0 0 0 1] * m
	for(int i = 0; i < 3; ++i)
		for(int j = 0; j < 4; ++j)
			m(i, j) *= s;
	return *this;
}

template<class T>
Transform<T> & Transform<T>::translate(const Array<T, 3> & v)
{
	// Equivalent to [1 0 0 v0 ; 0 1 0 v1 ; 0 0 1 v2 ; 0 0 0 1] * m
	m(0, 3) += v[0];
	m(1, 3) += v[1];
	m(2, 3) += v[2];

	return *this;
}

template<class T>
Transform<T> Transform<T>::combine_with(const Transform & t)
{
	Transform<T> ret;
	ret.m = m*t.m;
	return ret;
}


/* RigidParticle3D */
template<class T>
Box3D RigidParticle3D<T>::get_bounding_box() const
{
	return Box3D(
		position[0]-get_radius(), position[0]+get_radius(),
		position[1]-get_radius(), position[1]+get_radius(),
		position[2]-get_radius(), position[2]+get_radius()
	);
}

template<class T>
T RigidParticle3D<T>::compute_kinetic_energy() const
{
	const Array<T, 3> ang_vel_bf = orientation.apply_inv_rotation(ang_velocity);
	return 0.5*density*(normSqr(velocity) * shape->get_volume() +
			 	 	 	dot(shape->get_moment_of_inertia() * ang_vel_bf, ang_vel_bf));
}

template<class T>
void RigidParticle3D<T>::update()
{
	update_rotation_matrix();
	ang_momentum = compute_inertia_world_frame() * ang_velocity;
}

template<class T>
Matrix<T, 3> RigidParticle3D<T>::compute_inertia_world_frame() const
{
	return (rot_matrix.transpose() * shape->get_moment_of_inertia() * rot_matrix) * density * std::pow(scale, 5);
}

template<class T>
void RigidParticle3D<T>::get_node_info(pluint idx, ParticleNodeInfo<T> & ret) const
{
	const ParticleShapeNode<T> & node = shape->get_node(idx);
	ret.area = node.area * scale * scale;
	ret.normal = rot_matrix * node.normal;
	ret.pos_rel = rot_matrix * (scale * node.centroid);
	ret.pos = ret.pos_rel + position;
	ret.vel = velocity + crossProduct(ang_velocity, ret.pos_rel);
}

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
	/*const Quaternion<T> omega = Quaternion<T>(0, ang_velocity[0], ang_velocity[1], ang_velocity[2]);
	const Quaternion<T> omega_dot2 = Quaternion<T>(0, omega_dot[0], omega_dot[1], omega_dot[2]);
	const Quaternion<T> q_dot = 0.5 * omega;
	const Quaternion<T> q_dot_dot = 0.5 * (omega_dot2 * orientation + omega * q_dot);
	const T correction = 1 - 0.5*q_dot.norm_sqr() -
			std::sqrt(1 - q_dot.norm_sqr() - q_dot.dot(q_dot_dot) - 0.25*(q_dot_dot.norm_sqr() - util::sqr(q_dot.norm_sqr())));
	orientation = orientation + q_dot + 0.5*q_dot_dot - correction*orientation;*/

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

				/*const T projected_area = std::sqrt(shape->get_node(i).area * p2.shape->get_node(i).area) *
													dot(avg_n, get_normal(i)) *
													(-dot(avg_n, p2.get_normal(match)));
				Array<T, 3> df = -potential(dist) * projected_area * avg_n;*/

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

			if(std::abs(dist) <= potential.get_cutoff_distance()) {
				const Array<T, 3> df = -h * dot(h, normal) * node.area * potential(dist);
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
void RigidParticle3D<T>::write_to_stream_as_vtk(std::ostream & out) const
{
	Transform<T> transform = get_fixed_to_world_transform();

	out << "<?xml version=\"1.0\"?>\n";
	out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out << "  <UnstructuredGrid>\n";
	out << "    <Piece NumberOfPoints=\""<< shape->get_num_vertices() <<"\" NumberOfCells=\""<< shape->get_num_elements() <<"\">\n";
	out << "    <Points>\n";
	out << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
	Array<T, 3> p;
	const Array<T, 3> * vertices = shape->get_vertex_ptr();
	for(plint i = 0; i < shape->get_num_vertices(); ++i) {
		transform.apply(vertices[i], p);
		out << "        " << p[0] << " " << p[1] << " " << p[2] << std::endl;
	}
	out << "        </DataArray>\n";
	out << "    </Points>\n";
	out << "    <CellData>\n";
	out << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"contact\" Format=\"ascii\">\n";
	for(plint i = 0; i < collision_forces.size(); ++i) {
        out << std::setw(10) << std::setprecision(6)<< std::fixed
        	<< collision_forces[i][0] << " "<<collision_forces[i][1]<<" "<<collision_forces[i][2] << std::endl;
	}
	out << "        </DataArray>\n";
	out << "    </CellData>\n";
	out << "    <Cells>\n";
	out << "        <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n";
	const unsigned int * conn = shape->get_index_ptr();
	for(plint i = 0; i < shape->get_num_elements(); ++i) {
		out << "        " << conn[3*i] << " "
					<< conn[3*i+1] << " "
					<< conn[3*i+2] << std::endl;
	}

	out << "        </DataArray>\n";
	out << "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n";
	for(plint i=0; i < shape->get_num_elements(); i++)
		out << (i+1)*3 << " ";

	out << std::endl;
	out << "        </DataArray>\n";
	out << "        <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";

	for(int i=0; i < shape->get_num_elements(); i++)
		out << " 5";

	out << std::endl;

	out << "        </DataArray>\n";
	out << "      </Cells>\n";
	out << "    </Piece>\n";
	out << "  </UnstructuredGrid>\n";
	out << "</VTKFile>\n";
}

} /* namespace fsi */

} /* namespace plb */




#endif /* RIGIDPARTICLE_HH_ */
