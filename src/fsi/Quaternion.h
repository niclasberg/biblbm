/*
 * Quaternion.h
 *
 *  Created on: Dec 18, 2013
 *      Author: niber
 */

#ifndef QUATERNION_H_
#define QUATERNION_H_
#include <cmath>
#include "core/array.h"
#include "Matrix3.h"
#include <ostream>

namespace plb {

namespace fsi {

template<class T>
class Quaternion {
public:
	Quaternion()
	{
		_q[1] = _q[2] = _q[3] = (T) 0.0;
		_q[0] = (T) 1.0;
	}

	Quaternion(T angle, const Array<T, 3> & axis)
	{
		set_from_angle_axis(angle, axis);
	}

	void set_from_angle_axis(T angle, const Array<T, 3> & axis)
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

	T & operator[](unsigned int i)
	{
		return _q[i];
	}

	T operator[](unsigned int i) const
	{
		return _q[i];
	}

	T norm_sqr() const
	{
		return _q[0]*_q[0] + _q[1]*_q[1] + _q[2]*_q[2] + _q[3]*_q[3];
	}

	T norm() const
	{
		return std::sqrt(norm_sqr());
	}

	Quaternion normalize() const
	{
		T nrm = norm();
		return Quaternion<T>(_q[0]/nrm, _q[1]/nrm, _q[2]/nrm, _q[3]/nrm);
	}

	Quaternion conj() const
	{
		return Quaternion<T>(_q[0], -_q[1], -_q[2], -_q[3]);
	}

	T dot(const Quaternion<T> & q2) const {
		return _q[0]*q2[0] + _q[1]*q2[1] + _q[2]*q2[2] + _q[3]*q2[3];
	}

	void set_to_unity()
	{
		_q[1] = _q[2] = _q[3] = (T) 0.0;
		_q[0] = (T) 1.0;
	}

	void get_angle_axis(T & angle, Array<T, 3> & axis) const
	{
		T vec_norm = std::sqrt(_q[1]*_q[1] + _q[2]*_q[2] + _q[3]*_q[3]);
		angle = std::atan2(vec_norm, _q[0]);
		if(vec_norm < 1e-6)
			axis.set_to_zero();
		else {
			axis[0] = _q[1] / vec_norm;
			axis[1] = _q[2] / vec_norm;
			axis[2] = _q[3] / vec_norm;
		}
	}

	void apply_rotation(const Array<T, 3> & vec, Array<T, 3> & ret) const
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

	Array<T, 3> apply_rotation(const Array<T, 3> & vec) const
	{
		Array<T, 3> ret;
		apply_rotation(vec, ret);
		return ret;
	}

	void apply_inv_rotation(const Array<T, 3> & vec, Array<T, 3> & ret) const
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

	Array<T, 3> apply_inv_rotation(const Array<T, 3> & vec) const
	{
		Array<T, 3> ret;
		apply_inv_rotation(vec, ret);
		return ret;
	}

	void to_rot_matrix(Matrix<T, 3> & ret) const
	{
		ret(0, 0) = 1 - 2*(_q[2]*_q[2] + _q[3]*_q[3]);
		ret(0, 1) = 2*(_q[1]*_q[2] - _q[0]*_q[3]);
		ret(0, 2) = 2*(_q[0]*_q[2] + _q[1]*_q[3]);
		ret(1, 0) = 2*(_q[1]*_q[2] + _q[0]*_q[3]);
		ret(1, 1) = 1 - 2*(_q[1]*_q[1] + _q[3]*_q[3]);
		ret(1, 2) = 2*(_q[2]*_q[3] - _q[0]*_q[1]);
		ret(2, 0) = 2*(_q[1]*_q[3] - _q[0]*_q[2]);
		ret(2, 1) = 2*(_q[0]*_q[1] + _q[2]*_q[3]);
		ret(2, 2) = 1 - 2*(_q[1]*_q[1] + _q[2]*_q[2]);
	}

	Matrix<T, 3> to_rot_matrix() const
	{
		Matrix<T, 3> ret;
		to_rot_matrix(ret);
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

template<class T>
Quaternion<T> operator*(const Quaternion<T> & p, T s)
{
	return Quaternion<T>(s*p[0], s*p[1], s*p[2], s*p[3]);
}

template<class T>
Quaternion<T> operator*(T s, const Quaternion<T> & p)
{
	return Quaternion<T>(s*p[0], s*p[1], s*p[2], s*p[3]);
}

template<class T>
Quaternion<T> operator+(const Quaternion<T> & p, const Quaternion<T> & q)
{
	return Quaternion<T>(
				p[0]+q[0],
				p[1]+q[1],
				p[2]+q[2],
				p[3]+q[3]
		);
}

template<class T>
Quaternion<T> operator-(const Quaternion<T> & p, const Quaternion<T> & q)
{
	return Quaternion<T>(
				p[0]-q[0],
				p[1]-q[1],
				p[2]-q[2],
				p[3]-q[3]
		);
}

template<class T>
std::ostream & operator<<(std::ostream & ss, const Quaternion<T> & q)
{
	ss << "(" << q[0] << ", [" << q[1] << ", " << q[2] << ", " << q[3] << "])";
	return ss;
}

}

} /* namespace plb */

#endif /* QUATERNION_H_ */
