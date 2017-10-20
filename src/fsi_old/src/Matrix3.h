/*
 * Matrix3.h
 *
 *  Created on: Dec 18, 2013
 *      Author: niber
 */

#ifndef MATRIX3_H_
#define MATRIX3_H_
#include "Array3.h"


template<class T>
class Matrix3 {
public:
	Matrix3()
	{

	}

	T & operator()(unsigned int i, unsigned int j) {
		return _vals[i*3+j];
	}

	const T & operator()(unsigned int i, unsigned int j) const {
		return _vals[i*3+j];
	}

	Matrix3 transpose()
	{
		Matrix3<T> ret;
		ret(0, 0) = (*this)(0, 0);
		ret(0, 1) = (*this)(1, 0);
		ret(0, 2) = (*this)(2, 0);
		ret(1, 0) = (*this)(0, 1);
		ret(1, 1) = (*this)(1, 1);
		ret(1, 2) = (*this)(2, 1);
		ret(2, 0) = (*this)(0, 2);
		ret(2, 1) = (*this)(1, 2);
		ret(2, 2) = (*this)(2, 2);
		return ret;
	}

private:
	T _vals[9];
};

template<class T>
Array3<T> operator*(const Matrix3<T> & m, const Array3<T> & v)
{
	return Array3<T>(
				m(0, 0)*v[0] + m(0, 1)*v[1] + m(0, 2)*v[2],
				m(1, 0)*v[0] + m(1, 1)*v[1] + m(1, 2)*v[2],
				m(2, 0)*v[0] + m(2, 1)*v[1] + m(2, 2)*v[2]
			);
}

template<class T>
Array3<T> operator*(const Matrix3<T> & m, const Matrix3<T> & m2)
{
	Matrix3<T> ret;
	for(unsigned int i=0; i < 3; ++i) {
		for(unsigned int j=0; j < 3; ++j) {
			ret(i, j) = m(i, 0)*m2(0, j) + m(i, 1)*m2(1, j) + m(i, 2)*m2(2, j);
		}
	}
	return ret;
}

#endif /* MATRIX3_H_ */
