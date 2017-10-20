/*
 * Transform.h
 *
 *  Created on: Jun 2, 2015
 *      Author: niber
 */

#ifndef TRANSFORM_H_
#define TRANSFORM_H_
#include "Matrix3.h"
#include "Quaternion.h"

namespace plb {

namespace fsi {

template<class T>
struct Transform {
public:
	Transform();
	Array<T, 3> apply(const Array<T, 3> &) const;
	void apply(const Array<T, 3> &, Array<T, 3> &) const;
	Transform & rotate(const Matrix<T, 3> &);
	Transform & rotate(const Quaternion<T> &);
	Transform & scale(T);
	Transform combine_with(const Transform &);
	Transform & translate(const Array<T, 3> &);
private:
	Matrix<T, 4> m;
};

}

}




#endif /* TRANSFORM_H_ */
