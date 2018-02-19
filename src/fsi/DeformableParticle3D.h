/*
 * DeformableParticle3D.h
 *
 *  Created on: May 26, 2015
 *      Author: niber
 */

#ifndef DEFORMABLEPARTICLE3D_H_
#define DEFORMABLEPARTICLE3D_H_
#include "Quaternion.h"
#include "ParticleBase.h"

namespace plb {

namespace fsi {

template<class T>
class DeformableParticle3D : public ParticleBase3D<T> {
public:
	typedef typename ParticleBase3D<T>::vertex_const_iterator vertex_const_iterator;
	typedef typename ParticleBase3D<T>::vertex_iterator vertex_iterator;

	DeformableParticle3D(const ParticleShape<T> * shape);
	virtual ~DeformableParticle3D() { }

	virtual void move_vertices(Boundary<T> * boundary = 0);
	virtual void update();

	void compute_tensor_of_gyration(Matrix<T, 3> &) const;
	void compute_moment_of_intertia(Matrix<T, 3> &) const;
	void compute_orientation(Matrix<T, 3> &, Array<T, 3> &) const;

	// Transform vertices
	void transform_vertices(const Transform<T> &);
	void set_minor_axis_orientation(const Array<T, 3> &);
	void set_major_axis_orientation(const Array<T, 3> &);

	// IO
	virtual void write_lightweight(std::ostream &) const;
	virtual void pack(std::vector<char> &) const;
	virtual void unpack(char *&);
	virtual void unpack(std::istream &);

private:
	void rotate_to_align_vectors(const Array<T, 3> &, const Array<T, 3> &);
};

}

} /* namespace plb */

#endif /* DEFORMABLEPARTICLE3D_H_ */
