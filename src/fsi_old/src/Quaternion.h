/*
 * Quaternion.h
 *
 *  Created on: Dec 18, 2013
 *      Author: niber
 */

#ifndef QUATERNION_H_
#define QUATERNION_H_
#include <cmath>
#include "Array3.h"
#include "Matrix3.h"

template<class T>
class Quaternion {
public:
	Quaternion()
	{

	}

	Quaternion(T angle, const Array3<T> & axis)
	{
		T sin_a_half = std::sin(angle/2);
		_q[0] = std::cos(angle/2);
		_q[1] = sin_a_half*axis[0];
		_q[2] = sin_a_half*axis[1];
		_q[3] = sin_a_half*axis[2];
	}

	Quaternion(T q0, T q1, T q2, T q3)
	{
		_q[0] = q0;
		_q[1] = q1;
		_q[2] = q2;
		_q[3] = q3;
	}

	T norm_sqr() const
	{
		return _q[0]*_q[0] + _q[1]*_q[1] + _q[2]*_q[2] + _q[3]*_q[3];
	}

	T norm() const
	{
		return std::sqrt(norm_sqr());
	}

	Quaternion conj() const
	{
		return Quaternion<T>(_q[0], -_q[1], -_q[2], -_q[3]);
	}

	void apply_rotation(const Array3<T> & vec, Array3<T> & ret)
	{
		ret[0]= (1 - 2*(_q[2]*_q[2] + _q[3]*_q[3]))*vec[0] +
					  2*(_q[1]*_q[2] - _q[0]*_q[3])*vec[1] +
					  2*(_q[0]*_q[2] + _q[1]*_q[3])*vec[2];
		ret[1] = 2*(_q[1]*_q[2] + _q[0]*_q[3])*vec[0] +
					  (1 - 2*(_q[1]*_q[1] + _q[3]*_q[3]))*vec[1] +
					  2*(_q[2]*_q[3] - _q[0]*_q[1])*vec[2];
		ret[2] = 2*(_q[1]*_q[3] - _q[0]*_q[2])*vec[0] +
					  2*(_q[0]*_q[1] + _q[2]*_q[3])*vec[1] +
					  (1 - 2*(_q[1]*_q[1] + _q[2]*_q[2]))*vec[2];
	}

	Array3<T> apply_rotation(const Array3<T> & vec)
	{
		Array3<T> ret;
		apply_rotation(vec, ret);
		return ret;
	}

	void apply_inv_rotation(const Array3<T> & vec, Array3<T> & ret)
	{
		ret[0]= (1 - 2*(_q[2]*_q[2] + _q[3]*_q[3]))*vec[0] +
					  2*(_q[1]*_q[2] + _q[0]*_q[3])*vec[1] +
					  2*(_q[1]*_q[3] - _q[0]*_q[2])*vec[2];
		ret[1] = 2*(_q[1]*_q[2] - _q[0]*_q[3])*vec[0] +
					  (1 - 2*(_q[1]*_q[1] + _q[3]*_q[3]))*vec[1] +
					  2*(_q[0]*_q[1] + _q[2]*_q[3])*vec[2];
		ret[2] = 2*(_q[0]*_q[2] + _q[1]*_q[3])*vec[0] +
					  2*(_q[2]*_q[3] - _q[0]*_q[1])*vec[1] +
					  (1 - 2*(_q[1]*_q[1] + _q[2]*_q[2]))*vec[2];
	}

	Array3<T> apply_inv_rotation(const Array3<T> & vec)
	{
		Array3<T> ret;
		apply_inv_rotation(vec, ret);
		return ret;
	}

	Matrix3<T> to_rot_matrix() const
	{
		Matrix3<T> ret;
		ret(0, 0) = 1 - 2*(_q[2]*_q[2] + _q[3]*_q[3]);
		ret(0, 1) = 2*(_q[1]*_q[2] - _q[0]*_q[3]);
		ret(0, 2) = 2*(_q[0]*_q[2] + _q[1]*_q[3]);
		ret(1, 0) = 2*(_q[1]*_q[2] + _q[0]*_q[3]);
		ret(1, 1) = 1 - 2*(_q[1]*_q[1] + _q[3]*_q[3]);
		ret(1, 2) = 2*(_q[2]*_q[3] - _q[0]*_q[1]);
		ret(2, 0) = 2*(_q[1]*_q[3] - _q[0]*_q[2]);
		ret(2, 1) = 2*(_q[0]*_q[1] + _q[2]*_q[3]);
		ret(2, 2) = 1 - 2*(_q[1]*_q[1] + _q[2]*_q[2]);
		return ret;
	}

private:
	T _q[4];
};

template<class T>
Quaternion<T> operator*(const Quaternion<T> & p, const Quaternion<T> & q)
{
	return Quaternion<T>(
				p[0]*q[0]-(p[1]*q[1]+p[2]*q[2]+p[3]*q[3]),
				p[0]*q[1] + q[0]*p[1] + p[2]*q[3] - p[3]*q[2],
				p[0]*q[2] + q[0]*p[2] + p[3]*q[1] - p[1]*q[3],
				p[0]*q[3] + q[0]*p[3] + p[1]*q[2] - p[2]*q[1]
		);
}

#endif /* QUATERNION_H_ */
