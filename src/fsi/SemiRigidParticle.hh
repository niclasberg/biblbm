/*
 * SemiRigidParticle.hh
 *
 *  Created on: Sep 10, 2015
 *      Author: niber
 */

#ifndef SEMIRIGIDPARTICLE_HH_
#define SEMIRIGIDPARTICLE_HH_
#include "SemiRigidParticle.hh"

namespace plb {

namespace fsi {

/**** SemiRigidParticleParams ****/
template<class T>
	template<class BufferType>
void SemiRigidParticleParams<T>::pack(BufferType & buff) const
{
	utils::pack(buff, this->k_in_plane);
	utils::pack(buff, this->k_out_of_plane);
	utils::pack(buff, this->k_bend);
	utils::pack(buff, this->l0);
}

template<class T>
	template<class BufferType>
void SemiRigidParticleParams<T>::unpack(BufferType & buff)
{
	utils::unpack(buff, this->k_in_plane);
	utils::unpack(buff, this->k_out_of_plane);
	utils::unpack(buff, this->k_bend);
	utils::unpack(buff, this->l0);
}

/**** SemiRigidParticle3D ****/
template<class T>
plint SemiRigidParticle3D<T>::type_id = register_particle_type<T, SemiRigidParticle3D<T> >();

template<class T>
SemiRigidParticle3D<T>::SemiRigidParticle3D(const ParticleShape<T> * shape, const SemiRigidParticleParams<T> & params)
: DeformableParticle3D<T>(shape), params_(params)
{

}

template<class T>
SemiRigidParticle3D<T>::SemiRigidParticle3D(const ParticleShape<T> * shape)
: DeformableParticle3D<T>(shape), params_()
{

}

template<class T>
void SemiRigidParticle3D<T>::compute_forces()
{
	typedef typename ParticleShape<T>::link_const_iterator link_iterator;
	typedef typename ParticleShape<T>::triangle_const_iterator triangle_iterator;
	typedef typename ParticleShape<T>::triangle_pair_const_iterator triangle_pair_iterator;

	// Nodal forces (linear spring between the center of mass and each vertex)
	// The force on the center of mass position is evenly divided between the 
	// Nodes
	Array<T, 3> force_from_center_of_mass_displacement;
	force_from_center_of_mass_displacement.resetToZero();

	for(plint i = 0; i < this->count_nodes(); ++i) {
		Vertex<T> & v = this->get_node(i);
		const Array<T, 3> dx = v.pos - this->center_of_mass();
		const T l = norm(dx);
		const T l0 = norm(this->shape()->get_vertex(i) - this->shape()->get_center());

		const Array<T, 3> df = (-params().k_out_of_plane * (l - l0) / (l0*l)) * dx;
		v.force += df;
		force_from_center_of_mass_displacement -= df;
	}

	// Distribute the force on the center of mass to all nodes
	force_from_center_of_mass_displacement /= this->count_nodes();
	for(vertex_iterator it = this->begin(); it != this->end(); ++it)
		it->force += force_from_center_of_mass_displacement;

	// Link forces
	for(link_iterator it = this->shape()->links_begin();
			it != this->shape()->links_end(); ++it) {
		Vertex<T> & v0 = this->get_node(it->i0);
		Vertex<T> & v1 = this->get_node(it->i1);

		const Array<T, 3> dv = v0.pos - v1.pos;
		const T L = norm(dv);
		const T df = -params().k_in_plane * (L - it->length) / (L*it->length);
		v0.force += df * dv;
		v1.force -= df * dv;
	}
}

/******** Serialization *********/
template<class T>
void SemiRigidParticle3D<T>::pack(std::vector<char> & buff) const
{
	DeformableParticle3D<T>::pack(buff);
	params().pack(buff);
}

template<class T>
void SemiRigidParticle3D<T>::unpack(char *& buff)
{
	DeformableParticle3D<T>::unpack(buff);
	params().unpack(buff);
}

template<class T>
void SemiRigidParticle3D<T>::unpack(std::istream & buff)
{
	DeformableParticle3D<T>::unpack(buff);
	params().unpack(buff);
}

}

}


#endif /* SEMIRIGIDPARTICLE_HH_ */
