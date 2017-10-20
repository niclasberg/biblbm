/*
 * DeformableCapsule.hh
 *
 *  Created on: Dec 2, 2015
 *      Author: niber
 */

#ifndef DEFORMABLECAPSULE_HH_
#define DEFORMABLECAPSULE_HH_
#include "DeformableCapsule.h"
#include "ParticleFactory.h"
#include "utils.h"
#include "Potentials.h"
#include "TriangleUtils.h"

namespace plb {

namespace fsi {

/**** CapsuleParams ****/


template<class T>
	template<class BufferType>
void CapsuleParameters<T>::pack(BufferType & buff) const
{
	utils::pack(buff, this->C);
	utils::pack(buff, this->kb);
	utils::pack(buff, this->nu_s);
	utils::pack(buff, this->G);
}

template<class T>
	template<class BufferType>
void CapsuleParameters<T>::unpack(BufferType & buff)
{
	utils::unpack(buff, this->C);
	utils::unpack(buff, this->kb);
	utils::unpack(buff, this->nu_s);
	utils::unpack(buff, this->G);
}

/**** RBCParticle ****/
template<class T>
plint DeformableCapsuleParticle<T>::type_id = register_particle_type<T, DeformableCapsuleParticle<T> >();

template<class T>
DeformableCapsuleParticle<T>::DeformableCapsuleParticle(const ParticleShape<T> * shape, const CapsuleParameters<T> & params)
: DeformableParticle3D<T>(shape), params_(params), shear_model_(NeoHookean)
{

}

template<class T>
DeformableCapsuleParticle<T>::DeformableCapsuleParticle(const ParticleShape<T> * shape)
: DeformableParticle3D<T>(shape), params_(), shear_model_(NeoHookean)
{

}

template<class T>
void DeformableCapsuleParticle<T>::compute_forces()
{
	typedef typename ParticleShape<T>::triangle_const_iterator triangle_iterator;
	typedef typename ParticleShape<T>::triangle_quad_const_iterator triangle_quad_iterator;

	// Shear stress
	Array<Array<T, 3>, 3> grad_l, grad_lprime, grad_cos;

	for(triangle_iterator it = this->shape()->triangles_begin();
			it != this->shape()->triangles_end(); ++it) {
		Vertex<T> * v[3] = {&(this->get_node(it->i0)),
							&(this->get_node(it->i1)),
							&(this->get_node(it->i2))};

		// Edges in deformed and undeformed configuration
		const Array<T, 3> e0 = v[1]->pos - v[0]->pos;
		const Array<T, 3> e1 = v[2]->pos - v[0]->pos;
		const Array<T, 3> e0_0 = this->shape()->get_vertex(it->i1) - this->shape()->get_vertex(it->i0);
		const Array<T, 3> e1_0 = this->shape()->get_vertex(it->i2) - this->shape()->get_vertex(it->i0);

		// Edge lengths
		const T l = norm(e0);
		const T lprime = norm(e1);
		const T l0 = norm(e0_0);
		const T l0prime = norm(e1_0);
		const T e0_dot_e1 = dot(e0, e1);

		// Angles
		const T sin_phi = norm(crossProduct(e0, e1)) / (l*lprime);
		const T cos_phi = e0_dot_e1 / (l*lprime);
		const T sin_phi_0 = norm(crossProduct(e0_0, e1_0)) / (l0*l0prime);
		const T cos_phi_0 = dot(e0_0, e1_0) / (l0*l0prime);

		// Deformation gradient tensor components (a b; 0 c)
		const T a = l / l0;
		const T b = (lprime*cos_phi/l0prime - l*cos_phi_0/l0) / sin_phi_0;
		const T c = lprime * sin_phi / (lprime * sin_phi_0);

		// Green strain tensor
		const T D11 = (a*a - 1.) / 2.;
		const T D22 = (b*b + c*c - 1) / 2.;
		const T D12 = a*b/2.;

		// Invariants
		const T I1 = 2. * (D11 + D22);
		const T I2 = I1 + 4*(D11*D22 - D12*D12);

		// Gradient of edge 0 length
		grad_l[0] = e0 / (-l);
		grad_l[1] = e0 / l;
		grad_l[2].resetToZero();

		// Gradient of edge 1 length
		grad_lprime[0] = e1 / (-lprime);
		grad_lprime[1].resetToZero();
		grad_lprime[2] = e1 / lprime;

		// Gradient of cosine of angle between the edges
		grad_cos[0] = (-1. / (l*lprime) + e0_dot_e1 / (l*l*l*lprime)) * e0 +
					  (-1. / (l*lprime) + e0_dot_e1 / (l*lprime*lprime*lprime)) * e1;
		grad_cos[1] = (1. / (l*lprime)) * e1 - (e0_dot_e1 / (l*l*l*lprime)) * e0;
		grad_cos[2] = (1. / (l*lprime)) * e0 - (e0_dot_e1 / (l*lprime*lprime*lprime)) * e1;

		// Force = -grad(E) = -sum_{ij=1}^2 dW_s / dD_{ij} grad D_{ij}
		for(int i = 0; i < 3; ++i) {
			// Gradient of Green Strain tensor
			const Array<T, 3> grad_D11 = (l / (l0*l0)) * grad_l[i];
			const Array<T, 3> grad_D22 = (1. / util::sqr(sin_phi_0)) * (
							(lprime/util::sqr(l0prime) - l*cos_phi*cos_phi_0/(l0*l0prime)) * grad_lprime[i] +
							(l*util::sqr(cos_phi_0) / (l0*l0) - lprime*cos_phi*cos_phi_0/(l0prime*l0)) * grad_l[i] -
							(l*lprime*cos_phi_0 / (l0*l0prime)) * grad_cos[i]);
			const Array<T, 3> grad_D12 = (1. / (2.*l0prime*sin_phi_0)) * (
							(lprime*cos_phi/l0prime - 2.*l*cos_phi_0/l0) * grad_l[i] +
							(l*cos_phi/l0prime) * grad_lprime[i] +
							(l*lprime/l0prime) * grad_cos[i]);

			// Gradient of invariants
			const Array<T, 3> grad_I1 = 2. * (grad_D11 + grad_D22);
			const Array<T, 3> grad_I2 = (2. + 4.*D22)*grad_D11 + (2. + 4.*D11)*grad_D22 - 8.*D12*grad_D12;

			if(shear_model_ == Skalak) {
				// Skalak model: W_s = G/2 * (I1^2/2 + I1 - I2 + C*I2^2/2)
				v[i]->force += -params().G/2. * it->area * ((I1 + 1) * grad_I1 + (params().C * I2 - 1) * grad_I2);
			} else if(shear_model_ == NeoHookean) {
				// NeoHookean model: W_s = G/2 * (I1 - 1 + 1 / (I2 + 1))
				v[i]->force += -params().G/2. * it->area * (grad_I1 - (1. / util::sqr(1+I2)) * grad_I2 );
			} else {
				throw std::runtime_error("Unknown shear model!");
			}
		}
	}

	// Bending forces
	/*Array<Array<T, 3>, 6> grad_K11, grad_K12, grad_K22;

	for(triangle_quad_iterator it = this->shape()->triangle_quads_begin();
			it != this->shape()->triangle_quads_end(); ++it) {

		// Current vertex positions
		Vertex<T> * v[6] = {&(this->get_node(it->i0)),
						    &(this->get_node(it->i1)),
						    &(this->get_node(it->i2)),
						    &(this->get_node(it->i3)),
						    &(this->get_node(it->i4)),
						    &(this->get_node(it->i5))};

		// Vertex positions in undeformed configuration
		const Array<T, 3> & v00 = this->shape()->get_vertex(it->i0),
							v01 = this->shape()->get_vertex(it->i1),
							v02 = this->shape()->get_vertex(it->i2),
							v03 = this->shape()->get_vertex(it->i3),
							v04 = this->shape()->get_vertex(it->i4),
							v05 = this->shape()->get_vertex(it->i5);

		// Bending strain tensor
		const T K11_0 = norm(v[3]->pos - 2.*v[2]->pos + v[4]->pos);
		const T K11 = K11_0 - norm(v03 - 2.* v02 + v04);
		const T K12_0 = norm(v[0]->pos - v[1]->pos - v[2]->pos + v[3]->pos);
		const T K12 = K12_0 - norm(v00 - v01 - v02 + v03);
		const T K22_0 = norm(v[3]->pos - 2. * v[1]->pos + v[5]->pos);
		const T K22 = K22_0 - norm(v03 - 2. * v01 + v05);

		// Invariants
		const T J1 = K11 + K22;
		const T J2 = K11*K22 - K12*K12;

		// Gradient of Bending strain tensor
		grad_K11[0].resetToZero();
		grad_K11[1].resetToZero();
		grad_K11[2] = -2. * (v[3]->pos - 2.*v[2]->pos + v[4]->pos) / K11_0;
		grad_K11[3] = (v[3]->pos - 2.*v[2]->pos + v[4]->pos) / K11_0;
		grad_K11[4] = (v[3]->pos - 2.*v[2]->pos + v[4]->pos) / K11_0;
		grad_K11[5].resetToZero();

		grad_K22[0].resetToZero();
		grad_K22[1] = -2. * (v[3]->pos - 2. * v[1]->pos + v[5]->pos) / K22_0;
		grad_K22[2].resetToZero();
		grad_K22[3] = (v[3]->pos - 2. * v[1]->pos + v[5]->pos) / K22_0;
		grad_K22[4].resetToZero();
		grad_K22[5] = (v[3]->pos - 2. * v[1]->pos + v[5]->pos) / K22_0;

		grad_K12[0] = (v[0]->pos - v[1]->pos - v[2]->pos + v[3]->pos) / K12_0;
		grad_K12[1] = -(v[0]->pos - v[1]->pos - v[2]->pos + v[3]->pos) / K12_0;
		grad_K12[2] = -(v[0]->pos - v[1]->pos - v[2]->pos + v[3]->pos) / K12_0;
		grad_K12[3] = (v[0]->pos - v[1]->pos - v[2]->pos + v[3]->pos) / K12_0;
		grad_K12[4].resetToZero();
		grad_K12[5].resetToZero();


		for(int i = 0; i < 6; ++i) {
			v[i]->force += -params().kb/2. * (J1*(grad_K11[i] + grad_K22[i]) -
					(1-params().nu_s)*(K22 * grad_K11[i] - 2.*K12*grad_K12[i] + K11*grad_K22[i]));
		}
	}*/
}

template<class T>
void DeformableCapsuleParticle<T>::compute_wall_interaction_forces(const Boundary<T> & boundary)
{
	/*const MorsePotential<T> potential(2, 2, 1e-3);

	//bool ret = false;
	for(plint i = 0; i < this->count_nodes(); ++i) {
		Vertex<T> & node = this->get_node(i);

		T dist = boundary.distance_to_boundary(node.pos);

		if(dist <= potential.get_repulsion_distance()) {
			node.force += -potential(dist) * boundary.get_normal(node.pos);
			//ret = true;
		}
	}*/

	//return ret;
}

/******** Serialization *********/
template<class T>
void DeformableCapsuleParticle<T>::pack(std::vector<char> & buff) const
{
	DeformableParticle3D<T>::pack(buff);
	utils::pack(buff, shear_model_);
	params().pack(buff);
}

template<class T>
void DeformableCapsuleParticle<T>::unpack(char *& buff)
{
	DeformableParticle3D<T>::unpack(buff);
	utils::unpack(buff, shear_model_);
	params().unpack(buff);
}

template<class T>
void DeformableCapsuleParticle<T>::unpack(std::istream & buff)
{
	DeformableParticle3D<T>::unpack(buff);
	utils::unpack(buff, shear_model_);
	params().unpack(buff);
}

}

}



#endif /* DEFORMABLECAPSULE_HH_ */
