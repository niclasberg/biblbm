/*
 * DeformableCapsule.h
 *
 *  Created on: Nov 30, 2015
 *      Author: niber
 */

#ifndef DEFORMABLECAPSULE_H_
#define DEFORMABLECAPSULE_H_

#include "DeformableParticle3D.h"
#include "ParticleFactory.h"
#include <cmath>

namespace plb {

namespace fsi {

// Forward declarations
template<class T> class Boundary;

template<class T>
struct CapsuleParameters {
	template<class BufferType> void pack(BufferType &) const;
	template<class BufferType> void unpack(BufferType &);

	T G;
	T C;
	T nu_s;
	T kb;
};

template<class T>
class DeformableCapsuleParticle : public DeformableParticle3D<T> {
public:
	typedef typename DeformableParticle3D<T>::vertex_const_iterator vertex_const_iterator;
	typedef typename DeformableParticle3D<T>::vertex_iterator vertex_iterator;

	enum ShearModel {Skalak, NeoHookean};

	DeformableCapsuleParticle(const ParticleShape<T> *);
	DeformableCapsuleParticle(const ParticleShape<T> *, const CapsuleParameters<T> &);
	static DeformableCapsuleParticle * create(const ParticleShape<T> * shape) { return new DeformableCapsuleParticle(shape); }
	virtual ~DeformableCapsuleParticle() { }
	virtual DeformableCapsuleParticle * clone() const { return new DeformableCapsuleParticle(*this); };

	virtual plint get_type_id() const { return type_id; }

	virtual void compute_forces();
	void relax_nodes(T);
	void print_energies(T, std::ostream &) const;

	virtual void pack(std::vector<char> &) const;
	virtual void unpack(char *&);
	virtual void unpack(std::istream &);

	virtual bool should_voxelize() const { return true; }

	// Getters
	CapsuleParameters<T> & params() { return params_; }
	const CapsuleParameters<T> & params() const { return params_; }

private:
	ShearModel shear_model_;
	static plint type_id;
	CapsuleParameters<T> params_;
};

} /* namespace fsi */

} /* namespace plb */


#endif /* DEFORMABLECAPSULE_H_ */
