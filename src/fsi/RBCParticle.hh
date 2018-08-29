#ifndef RBCPARTICLE_HH_
#define RBCPARTICLE_HH_
#include "RBCParticle.h"
#include "ParticleFactory.h"
#include "utils.h"
#include "Potentials.h"
#include "TriangleUtils.h"

namespace plb {

namespace fsi {

/**** RBCParams ****/
template<class T> T RBCParameters<T>::x0(1. / 2.2);
template<class T> T RBCParameters<T>::x0_2(x0*x0);
template<class T> T RBCParameters<T>::x0_3(x0*x0*x0);
template<class T> T RBCParameters<T>::x0_4(x0*x0*x0*x0);
template<class T> T RBCParameters<T>::oneOverx0(2.2);

template<class T>
inline T RBCParameters<T>::lmax(const T & l0) const
{
	return l0 * oneOverx0;
}

template<class T>
void RBCParameters<T>::in_plane_ks_kp(T l0, T m, T & ks, T & kp) const
{
	// The relations for ks and kp are: (c.f. Fedosov 2010)
	// The energy should be minimal at the equilibrium length:
	//		0 = ks / 4 * (6x + (3x^2 - 2x^3)/(1-x)^2) - kp / l0^m
	// Shear modulus:
	//		4*l0 * G / sqrt(3) =  ks * ( x0/(2*(1-x0)^3) - 1/(4*(1-x0)^2) + 1/4) + kp * (m+1)/ l0^m

	static T ks_coeff1 = 0.25 * (6.*x0 + (3*x0*x0 - 2*x0*x0*x0)/((1-x0)*(1-x0)));
	static T ks_coeff2 = x0 / (2.*(1-x0)*(1-x0)*(1-x0)) + 0.25 * (1. - 1. / ((1-x0)*(1-x0)));

	const T kp_coeff1 = -1. / std::pow(l0, m);
	const T kp_coeff2 = (m+1) / std::pow(l0, m);
	const T G_lhs = (T)4. * l0 * shear_modulus / std::sqrt(3);

	ks = -kp_coeff1*G_lhs / (ks_coeff1*kp_coeff2 - ks_coeff2*kp_coeff1);
	kp =  ks_coeff1*G_lhs / (ks_coeff1*kp_coeff2 - ks_coeff2*kp_coeff1);
}

template<class T>
T RBCParameters<T>::cos_theta0() const
{
	return std::cos(this->theta0);
}

template<class T>
T RBCParameters<T>::sin_theta0() const
{
	return std::sin(this->theta0);
}

template<class T>
void RBCParameters<T>::in_plane_C(T & C) const
{
	static T factor = (4*x0*x0 - 9*x0 + 6) / ((1-x0)*(1-x0));
	static T A0 = std::sqrt(3) * L0*L0 / 4.;
	T ks;
	in_plane_ks(ks);
	C = std::sqrt(3) * A0*A0 * ks/ (4. * lmax(L0)) * factor;
}

template<class T>
void RBCParameters<T>::in_plane_ks(T & ks) const
{
	static T factor = (3. / (4.*(1-x0)*(1-x0)) - 0.75 + 4.*x0 + x0/(2*(1-x0)*(1-x0)*(1-x0))) / x0;
	ks = 4. * shear_modulus * lmax(L0) / (std::sqrt(3) * factor);
}

template<class T>
T RBCParameters<T>::K() const
{
	return 2 * shear_modulus + k_area_global + k_area_local;
}

template<class T>
T RBCParameters<T>::youngs_modulus() const
{
	return 4. * K() * shear_modulus / (K() + shear_modulus);
}

template<class T>
T RBCParameters<T>::poisson_ratio() const
{
	return (K() - shear_modulus) / (K() + shear_modulus);
}

template<class T>
	template<class BufferType>
void RBCParameters<T>::pack(BufferType & buff) const
{
	utils::pack(buff, this->L0);
	utils::pack(buff, this->shear_modulus);
	utils::pack(buff, this->k_area_global);
	utils::pack(buff, this->k_area_local);
	utils::pack(buff, this->k_bend);
	utils::pack(buff, this->k_volume);
	utils::pack(buff, this->theta0);
	utils::pack(buff, this->vol_desired);
}

template<class T>
	template<class BufferType>
void RBCParameters<T>::unpack(BufferType & buff)
{
	utils::unpack(buff, this->L0);
	utils::unpack(buff, this->shear_modulus);
	utils::unpack(buff, this->k_area_global);
	utils::unpack(buff, this->k_area_local);
	utils::unpack(buff, this->k_bend);
	utils::unpack(buff, this->k_volume);
	utils::unpack(buff, this->theta0);
	utils::unpack(buff, this->vol_desired);
}

/**** RBCParticle ****/
template<class T>
plint RBCParticle<T>::type_id = register_particle_type<T, RBCParticle<T> >();

template<class T>
RBCParticle<T>::RBCParticle(const ParticleShape<T> * shape, const RBCParameters<T> & params)
: DeformableParticle3D<T>(shape), params_(params)
{

}

template<class T>
RBCParticle<T>::RBCParticle(const ParticleShape<T> * shape)
: DeformableParticle3D<T>(shape), params_()
{

}

template<class T>
void RBCParticle<T>::compute_forces()
{
	typedef typename ParticleShape<T>::link_const_iterator link_iterator;
	typedef typename ParticleShape<T>::triangle_const_iterator triangle_iterator;
	typedef typename ParticleShape<T>::triangle_pair_const_iterator triangle_pair_iterator;
	static const T m = 2;

	// Contribution from link related forces
	T max_df_link = 0;
	{
		for(link_iterator it = this->shape()->links_begin(); it != this->shape()->links_end(); ++it) {
			// Compute spring properties
			T ks, kp;
			params().in_plane_ks_kp(it->length, m, ks, kp);

			// Compute node distance
			Vertex<T> & v0 = this->get_node(it->i0);
			Vertex<T> & v1 = this->get_node(it->i1);

			const Array<T, 3> dv = v0.pos - v1.pos;
			const T L = norm(dv);
			const T x = L / params().lmax(it->length);

			// Compute forces
			const T df = -0.25*ks*(6.*x + (3.*x*x - 2.*x*x*x)/((1-x)*(1-x)) ) + kp / std::pow(L, m);
			max_df_link = std::max(std::abs(df), max_df_link);
			v0.force += df * dv / L;
			v1.force -= df * dv / L;
		}
	}

	// Contribution from triangle level forces
	T local_area = T();
	T max_df_area = T();
	T max_df_vol = T();
	for(triangle_iterator it = this->shape()->triangles_begin(); it != this->shape()->triangles_end(); ++it) {
		Vertex<T> & v0 = this->get_node(it->i0);
		Vertex<T> & v1 = this->get_node(it->i1);
		Vertex<T> & v2 = this->get_node(it->i2);

		// Storage for gradients
		Array<T, 3> grad_v0, grad_v1, grad_v2;

		// Element area and area gradient
		tri::triangle_area_and_gradient(v0.pos, v1.pos, v2.pos, local_area, grad_v0, grad_v1, grad_v2);
		T df_area = 0.;

		// In-plane area compression contribution
		//T C;
		//params().in_plane_C(C);
		//df_area += C / local_area;

		// Global area constraint
		df_area += params().k_area_global * (this->shape()->get_area() - this->area()) / this->shape()->get_area();

		// Local area constraint
		df_area += params().k_area_local * (it->area - local_area) / it->area;

		// Find the nodal forces (chain-rule differentiation of the area terms)
		v0.force += df_area * grad_v0;
		v1.force += df_area * grad_v1;
		v2.force += df_area * grad_v2;

		// Volume constraint
		tri::grad_signed_volume(v0.pos, v1.pos, v2.pos, grad_v0, grad_v1, grad_v2);

		T df_vol = params().k_volume * (params().vol_desired - this->volume()) / params().vol_desired;
		v0.force += df_vol * grad_v0;
		v1.force += df_vol * grad_v1;
		v2.force += df_vol * grad_v2;

		max_df_area = std::max(std::abs(df_area), max_df_area);
		max_df_vol = std::max(std::abs(df_vol), max_df_vol);
	}

	//std::cout << "Link: " << max_df_link << ", Area: " << max_df_area << ", Volume: " << max_df_vol << std::endl;

	// Bending forces
	T cos_theta = T(), sin_theta = T();
	T cos_theta0 = params().cos_theta0(), sin_theta0 = params().sin_theta0();
	//Array<Array<T, 3>, 4> grad_cos_theta, grad_sin_theta;
	Array<Array<T, 3>, 4> grad_theta;
	for(triangle_pair_iterator it = this->shape()->triangle_pairs_begin();
			it != this->shape()->triangle_pairs_end(); ++it) {

		// The triangle pairs are stored such that the two triangles are
		// made up of (v0, v1, v2) and (v0, v2, v3)¸ respectively.
		Vertex<T> & v0 = this->get_node(it->i0);
		Vertex<T> & v1 = this->get_node(it->i1);
		Vertex<T> & v2 = this->get_node(it->i2);
		Vertex<T> & v3 = this->get_node(it->i3);

		// Evaluate cos and sine of the angle between normals
		tri::grad_angle_between_pair(v0.pos, v1.pos, v2.pos, v3.pos,
				cos_theta, sin_theta,
				grad_theta[0], grad_theta[1], grad_theta[2], grad_theta[3]);

		// The potential is U = kb(1 - cos(theta-theta0)) = kb*(1 - cos(theta)cos(theta0) - sin(theta)sin(theta0))
		// Hence, the force is F = -grad(U) = kb(cos(theta0)grad(cos(theta)) + sin(theta0)grad(sin(theta))) or
		// F = kb(-cos(theta0)sin(theta) + sin(theta0)cos(theta))grad(theta)
		//T prefactor = -params().k_bend * (cos_theta0 * sin_theta + sin_theta0 * cos_theta);
		T prefactor = params().k_bend * (-it->cos_theta0 * sin_theta + it->sin_theta0 * cos_theta);
		v0.force += prefactor * grad_theta[0];
		v1.force += prefactor * grad_theta[1];
		v2.force += prefactor * grad_theta[2];
		v3.force += prefactor * grad_theta[3];
	}
}

template<class T>
	template<class Stream>
void RBCParticle<T>::print_energies(Stream & out) const
{
	typedef typename ParticleShape<T>::link_const_iterator link_iterator;
	typedef typename ParticleShape<T>::triangle_const_iterator triangle_iterator;
	typedef typename ParticleShape<T>::triangle_pair_const_iterator triangle_pair_iterator;
	static const T m = 2;

	T E_in_plane = 0, E_bend = 0, E_vol = 0, E_area = 0;

	// Contribution from link related forces
	for(link_iterator it = this->shape()->links_begin();
			it != this->shape()->links_end(); ++it) {
		const Vertex<T> & v0 = this->get_node(it->i0);
		const Vertex<T> & v1 = this->get_node(it->i1);

		const T L = norm(v0.pos - v1.pos);
		const T Lmax = params().lmax(it->length);
		const T xi = L / Lmax;
		const T x0 = it->length / Lmax;

		T ks, kp;
		params().in_plane_ks_kp(it->length, m, ks, kp);

		// Evaluate energy (reference level set to zero at equilibrium)
		E_in_plane += ks*Lmax*((3*xi*xi - 2*xi*xi*xi) / (4.*(1 - xi)) - (3*x0*x0 - 2*x0*x0*x0) / (4.*(1 - x0))) + kp * (1. / ((m-1.) * std::pow(L, m-1)) - 1. / ((m-1.) * std::pow(it->length, m-1)));
	}

	// Contribution from triangle level forces
	for(triangle_iterator it = this->shape()->triangles_begin();
			it != this->shape()->triangles_end(); ++it) {
		// Element area
		T local_area = tri::triangle_area(this->get_node(it->i0).pos, this->get_node(it->i1).pos, this->get_node(it->i2).pos);

		// Area compression
		//E_in_plane += params().C / local_area;

		// Global area constraint
		E_area += params().k_area_global / (2. * this->shape()->get_area()) * util::sqr(this->area() - this->shape()->get_area());

		// Local area constraint
		E_area += params().k_area_local / (2. * it->area) * util::sqr(it->area - local_area);

		// Volume constraint
		E_vol += params().k_volume / (2 * params().vol_desired) * util::sqr(this->volume() - params_.vol_desired);
	}

	// Bending energy
	T sin_theta, cos_theta;
	T cos_theta0 = params().cos_theta0(), sin_theta0 = params().sin_theta0();
	for(triangle_pair_iterator it = this->shape()->triangle_pairs_begin();
			it != this->shape()->triangle_pairs_end(); ++it) {
	
		tri::cos_sin_of_angle_between_pair(
				this->get_node(it->i0).pos,
				this->get_node(it->i1).pos,
				this->get_node(it->i2).pos,
				this->get_node(it->i3).pos,
				cos_theta, sin_theta);
		//E_bend += params().k_bend * (1 - cos_theta*cos_theta0 - sin_theta*sin_theta0);
		E_bend += params().k_bend * (1 - cos_theta*it->cos_theta0 - sin_theta*it->sin_theta0);
	}

	out << "In plane: " << E_in_plane << ", Bending: " << E_bend << ", Volume constraint: " << E_vol << ", Area constraint: " << E_area << std::endl;
}


template<class T>
void RBCParticle<T>::relax_nodes(T damping)
{
	typedef typename ParticleShape<T>::triangle_const_iterator triangle_iterator;

	// Compute node masses
	std::vector<T> M(this->count_nodes(), T());
	for(triangle_iterator it = this->shape()->triangles_begin();
				it != this->shape()->triangles_end(); ++it) {
		T local_area = tri::triangle_area(this->get_node(it->i0).pos, this->get_node(it->i1).pos, this->get_node(it->i2).pos);
		M[it->i0] += local_area / 3.;
		M[it->i1] += local_area / 3.;
		M[it->i2] += local_area / 3.;
	}

	// Semi-implicit Euler integration
	plint ni = 0;
	for(vertex_iterator it = this->begin(); it != this->end(); ++it, ++ni) {
		it->vel += it->force / M[ni] - damping*it->vel;
		it->pos += it->vel;
	}

	this->update();
	this->reset_forces();
	this->compute_forces();
}

/******** Serialization *********/
template<class T>
void RBCParticle<T>::pack(std::vector<char> & buff) const
{
	DeformableParticle3D<T>::pack(buff);
	params().pack(buff);
}

template<class T>
void RBCParticle<T>::unpack(char *& buff)
{
	DeformableParticle3D<T>::unpack(buff);
	params().unpack(buff);
}

template<class T>
void RBCParticle<T>::unpack(std::istream & buff)
{
	DeformableParticle3D<T>::unpack(buff);
	params().unpack(buff);
}

}

}



#endif /* RBCPARTICLE_HH_ */
