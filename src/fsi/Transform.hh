/*
 * Transform.hh
 *
 *  Created on: Jun 2, 2015
 *      Author: niber
 */

#ifndef TRANSFORM_HH_
#define TRANSFORM_HH_
#include "Transform.h"

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
Transform<T> & Transform<T>::rotate(const Quaternion<T> & q)
{
	return rotate(q.to_rot_matrix());
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

}

}



#endif /* TRANSFORM_HH_ */
