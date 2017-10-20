/*
 * SemiRigidParticle.h
 *
 *  Created on: Sep 10, 2015
 *      Author: niber
 */

#ifndef SEMIRIGIDPARTICLE_H_
#define SEMIRIGIDPARTICLE_H_
#include "DeformableParticle3D.h"

namespace plb {

namespace fsi {

template<class T>
struct SemiRigidParticleParams {
	template<class BufferType> void pack(BufferType &) const;
	template<class BufferType> void unpack(BufferType &);

	T shear_modulus() const { return std::sqrt(3) * k_in_plane / (4 * l0); }

	// Average link length
	T l0;

	// In-plane shear
	T k_in_plane; 		// In-plane deformation stiffness

	// Out of plane properties
	T k_out_of_plane;

	// Bending parameters
	T k_bend;			// Bending resistance
};

template<class T>
class SemiRigidParticle3D : public DeformableParticle3D<T> {
public:
	typedef typename DeformableParticle3D<T>::vertex_const_iterator vertex_const_iterator;
	typedef typename DeformableParticle3D<T>::vertex_iterator vertex_iterator;

	SemiRigidParticle3D(const ParticleShape<T> *);
	SemiRigidParticle3D(const ParticleShape<T> *, const SemiRigidParticleParams<T> &);
	static SemiRigidParticle3D * create(const ParticleShape<T> * shape) { return new SemiRigidParticle3D(shape); }
	virtual ~SemiRigidParticle3D() { }
	virtual SemiRigidParticle3D * clone() const { return new SemiRigidParticle3D(*this); }

	virtual plint get_type_id() const { return type_id; }

	virtual void compute_forces();
	void print_energies(T, std::ostream &) const;

	virtual void pack(std::vector<char> &) const;
	virtual void unpack(char *&);
	virtual void unpack(std::istream &);

	// Getters
	SemiRigidParticleParams<T> & params() { return params_; }
	const SemiRigidParticleParams<T> & params() const { return params_; }

private:
	static plint type_id;
	SemiRigidParticleParams<T> params_;
};

} /* namespace fsi */

} /* namespace plb */


#endif /* SEMIRIGIDPARTICLE_H_ */
