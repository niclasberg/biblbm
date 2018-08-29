/*
 * RBCParticle.h
 *
 *  Created on: Jun 1, 2015
 *      Author: niber
 */

#ifndef RBCPARTICLE_H_
#define RBCPARTICLE_H_
#include "DeformableParticle3D.h"
#include "ParticleFactory.h"
#include <cmath>

namespace plb {

namespace fsi {

// Forward declarations
template<class T> class Boundary;

template<class T>
struct RBCParameters {
	T G() const { return shear_modulus; }
	T K() const;
	T youngs_modulus() const;
	T poisson_ratio() const;

	T lmax(const T &) const;
	void in_plane_C(T & C) const;
	void in_plane_ks(T & ks) const;
	void in_plane_ks_kp(T l0, T m, T & ks, T & kp) const;

	T cos_theta0() const;
	T sin_theta0() const;

	template<class BufferType> void pack(BufferType &) const;
	template<class BufferType> void unpack(BufferType &);

	// Non-dimensional parameters
	T shear_modulus;
	T L0;

	// Area constraint
	T k_area_global;
	T k_area_local;		// Local area constraint stiffness

	// Volume constraint
	T k_volume;			// Global volume constraint stiffness
	T vol_desired;		// Desired volume of the cell

	// Bending parameters
	T k_bend;			// Bending resistance
	T theta0;			// Spontaneous angle

private:
	static T sqrt3;
	static T x0, x0_2, x0_3, x0_4, oneOverx0;
};

template<class T>
class RBCParticle : public DeformableParticle3D<T> {
public:
	typedef typename DeformableParticle3D<T>::vertex_const_iterator vertex_const_iterator;
	typedef typename DeformableParticle3D<T>::vertex_iterator vertex_iterator;

	RBCParticle(const ParticleShape<T> *);
	RBCParticle(const ParticleShape<T> *, const RBCParameters<T> &);
	static RBCParticle * create(const ParticleShape<T> * shape) { return new RBCParticle(shape); }
	virtual ~RBCParticle() { }
	virtual RBCParticle * clone() const { return new RBCParticle(*this); };

	virtual plint get_type_id() const { return type_id; }
	virtual void compute_forces();
	void relax_nodes(T);
	template<class Stream> void print_energies(Stream &) const;

	virtual void pack(std::vector<char> &) const;
	virtual void unpack(char *&);
	virtual void unpack(std::istream &);

	// Getters
	RBCParameters<T> & params() { return params_; }
	const RBCParameters<T> & params() const { return params_; }

private:
	static plint type_id;
	RBCParameters<T> params_;
};

}

}



#endif /* RBCPARTICLE_H_ */
