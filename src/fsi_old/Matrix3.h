/*
 * Matrix3.h
 *
 *  Created on: Dec 18, 2013
 *      Author: niber
 */

#ifndef MATRIX3_H_
#define MATRIX3_H_
#include "core/array.h"
#include <algorithm>

namespace plb {

namespace fsi {

// Arbitrary square matrices (not implemented)
template<class T, int dim>
class Matrix {
public:

private:
	T _vals[dim*dim];
};

template<class T, int dim>
Matrix<T, dim> operator*(const Matrix<T, dim> & m, const Matrix<T, dim> & m2)
{
	Matrix<T, dim> ret;
	for(unsigned int i=0; i < dim; ++i) {
		for(unsigned int j=0; j < dim; ++j) {
			ret(i, j) = 0;
			for(unsigned int k = 0; k < dim; ++k)
				ret(i, j) += m(i, k)*m2(k, j);
		}
	}
	return ret;
}

// 3x3 matrix
template<class T>
class Matrix<T, 3> {
public:
	Matrix()
	{
		_vals[0] = T();
		_vals[1] = T();
		_vals[2] = T();
		_vals[3] = T();
		_vals[4] = T();
		_vals[5] = T();
		_vals[6] = T();
		_vals[7] = T();
		_vals[8] = T();
	}

	Matrix(const Matrix & rhs)
	{
		std::copy(rhs._vals, rhs._vals+9, _vals);
	}

	Matrix & operator=(const Matrix<T, 3> & rhs)
	{
		std::copy(rhs._vals, rhs._vals+9, _vals);
		return *this;
	}

	T & operator[](unsigned int i)
	{
		return _vals[i];
	}

	T operator[](unsigned int i) const
	{
		return _vals[i];
	}

	T & operator()(unsigned int i, unsigned int j) {
		return _vals[i*3+j];
	}

	const T & operator()(unsigned int i, unsigned int j) const {
		return _vals[i*3+j];
	}

	Matrix transpose() const
	{
		Matrix<T, 3> ret;
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

	Array<T, 3> transpose_multiply(const Array<T, 3> & rhs) const
	{
		return Array<T, 3>(
							rhs[0]*_vals[0] + rhs[1]*_vals[3] + rhs[2]*_vals[6],
							rhs[0]*_vals[1] + rhs[1]*_vals[4] + rhs[2]*_vals[7],
							rhs[0]*_vals[2] + rhs[1]*_vals[5] + rhs[2]*_vals[8]
						);
	}

	Array<T, 3> solve(const Array<T, 3> & rhs) const
	{
		return Array<T, 3>(
				rhs[2]*(_vals[1]*_vals[5] - _vals[2]*_vals[4]) +
					rhs[1]*(_vals[2]*_vals[7] - _vals[1]*_vals[8]) +
					rhs[0]*(_vals[4]*_vals[8] - _vals[5]*_vals[7]),
				rhs[2]*(_vals[2]*_vals[3] - _vals[0]*_vals[5]) +
					rhs[1]*(_vals[0]*_vals[8] - _vals[2]*_vals[6]) +
					rhs[0]*(_vals[5]*_vals[6] - _vals[3]*_vals[8]),
				rhs[2]*(_vals[0]*_vals[4] - _vals[1]*_vals[3])  +
					rhs[1]*(_vals[1]*_vals[6] - _vals[0]*_vals[7]) +
					rhs[0]*(_vals[3]*_vals[7] - _vals[4]*_vals[6])
			) / (_vals[0]*_vals[4]*_vals[8] - _vals[0]*_vals[5]*_vals[7] -
				 _vals[1]*_vals[3]*_vals[8] + _vals[1]*_vals[5]*_vals[6] +
				 _vals[2]*_vals[3]*_vals[7] - _vals[2]*_vals[4]*_vals[6]);
	}

	T * data()
	{
		return _vals;
	}

	const T * data() const
	{
		return _vals;
	}

private:
	T _vals[9];
};

template<class T>
Array<T, 3> operator*(const Matrix<T, 3> & m, const Array<T, 3> & v)
{
	return Array<T, 3>(
				m(0, 0)*v[0] + m(0, 1)*v[1] + m(0, 2)*v[2],
				m(1, 0)*v[0] + m(1, 1)*v[1] + m(1, 2)*v[2],
				m(2, 0)*v[0] + m(2, 1)*v[1] + m(2, 2)*v[2]
			);
}

template<class T>
Matrix<T, 3> operator*(const Matrix<T, 3> & m, const Matrix<T, 3> & m2)
{
	Matrix<T, 3> ret;
	for(unsigned int i=0; i < 3; ++i) {
		for(unsigned int j=0; j < 3; ++j) {
			ret(i, j) = m(i, 0)*m2(0, j) + m(i, 1)*m2(1, j) + m(i, 2)*m2(2, j);
		}
	}
	return ret;
}

template<class T>
Matrix<T, 3> operator*(const Matrix<T, 3> & m, const T & s)
{
	Matrix<T, 3> ret;
	for(unsigned int i=0; i < 9; ++i)
		ret[i] = s * m[i];
	return ret;
}

template<class T>
Matrix<T, 3> operator*(const T & s, const Matrix<T, 3> & m)
{
	Matrix<T, 3> ret;
	for(unsigned int i=0; i < 9; ++i)
		ret[i] = s * m[i];
	return ret;
}


// 4x4 matrix
template<class T>
class Matrix<T, 4> {
public:
	Matrix()
	{
		for(unsigned int i = 0; i< 16; ++i)
			_vals[i] = T();
	}

	Matrix(const Matrix & rhs)
	{
		std::copy(rhs._vals, rhs._vals+16, _vals);
	}

	Matrix & operator=(const Matrix & rhs)
	{
		std::copy(rhs._vals, rhs._vals+16, _vals);
		return *this;
	}

	T & operator[](unsigned int i)
	{
		return _vals[i];
	}

	T operator[](unsigned int i) const
	{
		return _vals[i];
	}

	T & operator()(unsigned int i, unsigned int j) {
		return _vals[i*4+j];
	}

	const T & operator()(unsigned int i, unsigned int j) const {
		return _vals[i*4+j];
	}

	Matrix transpose()
	{
		Matrix<T, 3> ret;
		ret(0, 0) = (*this)(0, 0);
		ret(0, 1) = (*this)(1, 0);
		ret(0, 2) = (*this)(2, 0);
		ret(0, 3) = (*this)(3, 0);
		ret(1, 0) = (*this)(0, 1);
		ret(1, 1) = (*this)(1, 1);
		ret(1, 2) = (*this)(2, 1);
		ret(1, 3) = (*this)(3, 1);
		ret(2, 0) = (*this)(0, 2);
		ret(2, 1) = (*this)(1, 2);
		ret(2, 2) = (*this)(2, 2);
		ret(2, 3) = (*this)(3, 2);
		ret(3, 0) = (*this)(0, 3);
		ret(3, 1) = (*this)(1, 3);
		ret(3, 2) = (*this)(2, 3);
		ret(3, 3) = (*this)(3, 3);
		return ret;
	}

	T * data()
	{
		return _vals;
	}

	const T * data() const
	{
		return _vals;
	}

private:
	T _vals[16];
};

template<class T>
Matrix<T, 4> operator*(const Matrix<T, 4> & m, const Matrix<T, 4> & m2)
{
	Matrix<T, 4> ret;
	for(unsigned int i=0; i < 4; ++i) {
		for(unsigned int j=0; j < 4; ++j) {
			ret(i, j) = m(i, 0)*m2(0, j) + m(i, 1)*m2(1, j) + m(i, 2)*m2(2, j) + m(i, 3)*m2(3, j);
		}
	}
	return ret;
}

}

} /* namespace plb */

#endif /* MATRIX3_H_ */
