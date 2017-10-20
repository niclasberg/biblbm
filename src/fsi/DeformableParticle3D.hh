/*
 * DeformableParticle3D.hh
 *
 *  Created on: May 26, 2015
 *      Author: niber
 */

#ifndef DEFORMABLEPARTICLE3D_HH_
#define DEFORMABLEPARTICLE3D_HH_

#include "DeformableParticle3D.h"
#include "linalg.h"

namespace plb {

namespace fsi {

template<class T>
DeformableParticle3D<T>::DeformableParticle3D(const ParticleShape<T> * shape)
: ParticleBase3D<T>(shape)
{
	this->update();
}

template<class T>
void DeformableParticle3D<T>::move_vertices()
{
	// Integrate with Explicit Euler
	for(vertex_iterator it = this->begin(); it != this->end(); ++it)
		it->pos += it->vel;
	update();
}

template<class T>
void DeformableParticle3D<T>::update()
{
	ParticleBase3D<T>::update();

	typedef typename ParticleShape<T>::triangle_const_iterator triangle_iterator;

	// Calculate total volume and area
	this->area_ = 0;
	this->volume_ = 0;
	this->center_of_mass_.resetToZero();

	Array<T, 3> normal, centroid;
	T dA;
	for(triangle_iterator it = this->shape()->triangles_begin();
			it != this->shape()->triangles_end(); ++it) {
		Array<T, 3> v0 = this->get_node(it->i0).pos;
		Array<T, 3> v1 = this->get_node(it->i1).pos;
		Array<T, 3> v2 = this->get_node(it->i2).pos;

		tri::centroid(v0, v1, v2, centroid);
		tri::normal_and_area(v0, v1, v2, normal, dA);

		this->area_ += dA;
		this->volume_ += dot(centroid, normal) * dA / 3.0;
		this->center_of_mass_ += centroid * dA;
	}

	this->center_of_mass_ /= this->area_;
}

template<class T>
void DeformableParticle3D<T>::set_minor_axis_orientation(const Array<T, 3> & pref_dir)
{
	// Compute orientation of the minor axis
	Matrix<T, 3> Q;
	Array<T, 3> eigs;
	this->compute_orientation(Q, eigs);

	// Smallest eigenvalue (i.e. minor axis)
	Array<T, 3> a_eigs(std::abs(eigs[0]), std::abs(eigs[1]), std::abs(eigs[2]));
	int k = (a_eigs[0] < a_eigs[1] && a_eigs[0] < a_eigs[2]) ? 0 : (a_eigs[1] < a_eigs[2]) ? 1 : 2;
	Array<T, 3> dir(Q(0, k), Q(1, k), Q(2, k));
	dir /= norm(dir);

	rotate_to_align_vectors(dir, pref_dir);
}

template<class T>
void DeformableParticle3D<T>::compute_orientation(Matrix<T, 3> & Q, Array<T, 3> & eig_values) const
{
	// Compute the tensor of gyration
	Matrix<T, 3> A;
	compute_tensor_of_gyration(A);
	linalg::diagonalize(A, Q, eig_values);
}

template<class T>
void DeformableParticle3D<T>::compute_tensor_of_gyration(Matrix<T, 3> & A) const
{
	A.reset_to_zero();
	for(vertex_const_iterator it = this->begin(); it != this->end(); ++it) {
		for(plint i = 0; i < 3; ++i)
			for(plint j = i; j < 3; ++j)
				A(i, j) += (it->pos[i] - this->center_of_mass()[i])*(it->pos[j] - this->center_of_mass()[j]);
	}
	A(1, 0) = A(0, 1);
	A(2, 0) = A(0, 2);
	A(2, 1) = A(1, 2);

	for(plint i = 0; i < 9; ++i)
		A[i] /= this->count_nodes();
}

template<class T>
void DeformableParticle3D<T>::rotate_to_align_vectors(const Array<T, 3> & dir, const Array<T, 3> & pref_dir)
{
	Array<T, 3> rot_axis = crossProduct(dir, pref_dir);
	T angle = std::asin(norm(rot_axis) / (norm(dir) * norm(pref_dir)));

	rot_axis /= norm(rot_axis);

	// Make sure that the sign of the angle is correct
	if(dot(pref_dir, dir) < 0)
		rot_axis = -rot_axis;

	Transform<T> transform;
	transform.translate(-this->center_of_mass())
			 .rotate(Quaternion<T>(angle, rot_axis))
			 .translate(this->center_of_mass());
	this->transform_vertices(transform);
}

template<class T>
void DeformableParticle3D<T>::transform_vertices(const Transform<T> & transform)
{
	for(vertex_iterator it = this->begin(); it != this->end(); ++it) {
		it->pos = transform.apply(it->pos);
	}
	this->update();
}

template<class T>
void DeformableParticle3D<T>::set_major_axis_orientation(const Array<T, 3> & pref_dir)
{
	// Compute orientation of the minor axis
	Matrix<T, 3> Q;
	Array<T, 3> eigs;
	compute_orientation(Q, eigs);

	// Largest eigenvalue (i.e. major axis)
	Array<T, 3> a_eigs(std::abs(eigs[0]), std::abs(eigs[1]), std::abs(eigs[2]));
	int k = (a_eigs[0] > a_eigs[1] && a_eigs[0] > a_eigs[2]) ? 0 : (a_eigs[1] > a_eigs[2]) ? 1 : 2;
	Array<T, 3> dir(Q(0, k), Q(1, k), Q(2, k));
	dir /= norm(dir);

	rotate_to_align_vectors(dir, pref_dir);
}

template<class T>
void DeformableParticle3D<T>::compute_moment_of_intertia(Matrix<T, 3> & moment_of_inertia) const
{
	typedef typename ParticleShape<T>::triangle_const_iterator triangle_iterator;

	const T oneDiv6 = (T)(1.0/6.0);
	const T oneDiv24 = (T)(1.0/24.0);
	const T oneDiv60 = (T)(1.0/60.0);
	const T oneDiv120 = (T)(1.0/120.0);

	// order:  1, x, y, z, x^2, y^2, z^2, xy, yz, zx
	T integral[10] = { (T)0.0, (T)0.0, (T)0.0, (T)0.0,
		(T)0.0, (T)0.0, (T)0.0, (T)0.0, (T)0.0, (T)0.0 };

	for(triangle_iterator it = this->shape()->triangles_begin(); it != this->shape()->triangles_end(); ++it) {
		const Array<T, 3> & v0 = this->get_node(it->i0).pos;
		const Array<T, 3> & v1 = this->get_node(it->i1).pos;
		const Array<T, 3> & v2 = this->get_node(it->i2).pos;

		// Get cross product of edges and normal vector.
		Array<T, 3> V1mV0 = v1 - v0;
		Array<T, 3> V2mV0 = v2 - v0;
		Array<T, 3> N = crossProduct(V1mV0, V2mV0);

		// Compute integral terms.
		T tmp0, tmp1, tmp2;
		T f1x, f2x, f3x, g0x, g1x, g2x;
		tmp0 = v0[0] + v1[0];
		f1x = tmp0 + v2[0];
		tmp1 = v0[0]*v0[0];
		tmp2 = tmp1 + v1[0]*tmp0;
		f2x = tmp2 + v2[0]*f1x;
		f3x = v0[0]*tmp1 + v1[0]*tmp2 + v2[0]*f2x;
		g0x = f2x + v0[0]*(f1x + v0[0]);
		g1x = f2x + v1[0]*(f1x + v1[0]);
		g2x = f2x + v2[0]*(f1x + v2[0]);

		T f1y, f2y, f3y, g0y, g1y, g2y;
		tmp0 = v0[1] + v1[1];
		f1y = tmp0 + v2[1];
		tmp1 = v0[1]*v0[1];
		tmp2 = tmp1 + v1[1]*tmp0;
		f2y = tmp2 + v2[1]*f1y;
		f3y = v0[1]*tmp1 + v1[1]*tmp2 + v2[1]*f2y;
		g0y = f2y + v0[1]*(f1y + v0[1]);
		g1y = f2y + v1[1]*(f1y + v1[1]);
		g2y = f2y + v2[1]*(f1y + v2[1]);

		T f1z, f2z, f3z, g0z, g1z, g2z;
		tmp0 = v0[2] + v1[2];
		f1z = tmp0 + v2[2];
		tmp1 = v0[2]*v0[2];
		tmp2 = tmp1 + v1[2]*tmp0;
		f2z = tmp2 + v2[2]*f1z;
		f3z = v0[2]*tmp1 + v1[2]*tmp2 + v2[2]*f2z;
		g0z = f2z + v0[2]*(f1z + v0[2]);
		g1z = f2z + v1[2]*(f1z + v1[2]);
		g2z = f2z + v2[2]*(f1z + v2[2]);

		// Update integrals.
		integral[0] += N[0]*f1x;
		integral[1] += N[0]*f2x;
		integral[2] += N[1]*f2y;
		integral[3] += N[2]*f2z;
		integral[4] += N[0]*f3x;
		integral[5] += N[1]*f3y;
		integral[6] += N[2]*f3z;
		integral[7] += N[0]*(v0[1]*g0x + v1[1]*g1x + v2[1]*g2x);
		integral[8] += N[1]*(v0[2]*g0y + v1[2]*g1y + v2[2]*g2y);
		integral[9] += N[2]*(v0[0]*g0z + v1[0]*g1z + v2[0]*g2z);
	}

	integral[0] *= oneDiv6;
	integral[1] *= oneDiv24;
	integral[2] *= oneDiv24;
	integral[3] *= oneDiv24;
	integral[4] *= oneDiv60;
	integral[5] *= oneDiv60;
	integral[6] *= oneDiv60;
	integral[7] *= oneDiv120;
	integral[8] *= oneDiv120;
	integral[9] *= oneDiv120;

	// center of mass
	T volume = integral[0];
	Array<T, 3> center(integral[1]/volume, integral[2]/volume, integral[3]/volume);

	// inertia relative to center of mass
	moment_of_inertia(0, 0) = integral[5] + integral[6] -
			volume*(center[1]*center[1] + center[2]*center[2]);
	moment_of_inertia(0, 1) = -integral[7] + volume*center[0]*center[1];
	moment_of_inertia(0, 2) = -integral[9] + volume*center[2]*center[0];
	moment_of_inertia(1, 0) = moment_of_inertia(0, 1);
	moment_of_inertia(1, 1) = integral[4] + integral[6] -
			volume*(center[2]*center[2] + center[0]*center[0]);
	moment_of_inertia(1, 2) = -integral[8] + volume*center[1]*center[2];
	moment_of_inertia(2, 0) = moment_of_inertia(0, 2);
	moment_of_inertia(2, 1) = moment_of_inertia(1, 2);
	moment_of_inertia(2, 2) = integral[4] + integral[5] -
			volume*(center[0]*center[0] + center[1]*center[1]);
}


/*********** IO ************/
template<class T>
void DeformableParticle3D<T>::pack(std::vector<char> & buff) const
{
	ParticleBase3D<T>::pack(buff);

	// Pack vertices
	for(vertex_const_iterator it = this->begin(); it != this->end(); ++it) {
		utils::pack(buff, it->pos);
		utils::pack(buff, it->vel);
		utils::pack(buff, it->force);
	}
}

template<class T>
void DeformableParticle3D<T>::unpack(char *& buff)
{
	ParticleBase3D<T>::unpack(buff);

	// Unpack vertices
	for(vertex_iterator it = this->begin(); it != this->end(); ++it) {
		utils::unpack(buff, it->pos);
		utils::unpack(buff, it->vel);
		utils::unpack(buff, it->force);
	}
	this->update();
}

template<class T>
void DeformableParticle3D<T>::unpack(std::istream & in)
{
	ParticleBase3D<T>::unpack(in);

	// Unpack vertices
	for(vertex_iterator it = this->begin(); it != this->end(); ++it) {
		utils::unpack(in, it->pos);
		utils::unpack(in, it->vel);
		utils::unpack(in, it->force);
	}
	this->update();
}

template<class T>
void DeformableParticle3D<T>::write_lightweight(std::ostream & out) const
{
	out << this->shape()->get_tag() << ":" << this->get_id();

	// Print center of mass
	for(plint i = 0; i < 3; ++i)
		out << ":" << this->center_of_mass()[i];

	// Print tensor of gyration
	Matrix<T, 3> gyration_tensor;
	this->compute_tensor_of_gyration(gyration_tensor);
	for(plint i = 0; i < 3; ++i)
		for(plint j = i; j < 3; ++j)
			out << ":" << gyration_tensor(i, j);
	out << ":" << this->volume();
}



} /* namespace fsi */

} /* namespace plb */

#endif /* DEFORMABLEPARTICLE3D_HH_ */
