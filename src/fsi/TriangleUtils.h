#ifndef TRIANGLEUTILS_H_
#define TRIANGLEUTILS_H_
#include "core/array.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include <cmath>

namespace plb {

namespace fsi {

namespace tri {

template<class T>
inline void centroid(const Array<T, 3> & v0, const Array<T, 3> & v1, const Array<T, 3> & v2, Array<T, 3> & res)
{
	res = (v0 + v1 + v2) / (T) 3;
}

template<class T>
inline Array<T, 3> centroid(const Array<T, 3> & v0, const Array<T, 3> & v1, const Array<T, 3> & v2)
{
	Array<T, 3> res;
	centroid(v0, v1, v2, res);
	return res;
}

template<class T>
inline void normal_and_area(const Array<T, 3> & v0, const Array<T, 3> & v1, const Array<T, 3> & v2, Array<T, 3> & res, T & area)
{
	crossProduct(v1-v0, v2-v0, res);
	area = norm(res);
	res /= area;
	area *= 0.5;
}

template<class T>
inline void normal(const Array<T, 3> & v0, const Array<T, 3> & v1, const Array<T, 3> & v2, Array<T, 3> & res)
{
	crossProduct(v1-v0, v2-v0, res);
	res /= norm(res);
}

template<class T>
inline Array<T, 3> normal(const Array<T, 3> & v0, const Array<T, 3> & v1, const Array<T, 3> & v2)
{
	Array<T, 3> res;
	normal(v0, v1, v2, res);
	return res;
}

template<class T>
inline T triangle_area(const Array<T, 3> & v0, const Array<T, 3> & v1, const Array<T, 3> & v2)
{
	return 0.5*norm(crossProduct(v1-v0, v2-v0));
}

template<class T>
inline void triangle_area_and_gradient(const Array<T, 3> & v0, const Array<T, 3> & v1, const Array<T, 3> & v2,
		T & area, Array<T, 3> & g0, Array<T, 3> & g1, Array<T, 3> & g2)
{
	// Precompute some differences that are commonly occurring
	const T v01v11 = v0[1] - v1[1];
	const T v02v22 = v0[2] - v2[2];
	const T v02v12 = v0[2] - v1[2];
	const T v01v21 = v0[1] - v2[1];
	const T v00v20 = v0[0] - v2[0];
	const T v00v10 = v0[0] - v1[0];

	// Evaluate components (a, b, c) = crossProduct(v1-v0, v2-v0):
	const T a = v01v11*v02v22 - v02v12*v01v21;
	const T b = v02v12*v00v20 - (v0[0] - v1[0])*v02v22;
	const T c = v00v10*v01v21 - v01v11*v00v20;

	// The area is given by 0.5 * norm(crossProduct(v1-v0, v2-v0))
	// or, sqrt(a^2 + b^2 + c^2) / 2
	area = 0.5*std::sqrt(a*a + b*b + c*c);

	// Will use that grad(A) = grad(A^2) / (2*A)
	// Now, grad(A^2) = grad(a^2 + b^2 + c^2) / 4 = (a grad(a) + b grad(b) + c grad(c)) / 2
	// Thus grad(A) = (a grad(a) + b grad(b) + c grad(c)) / (4*A)
	const T denom = (T) 1. / ((T)4. * area);

	//	grad_v0(a) = (0, v1[2] - v2[2], v2[1] - v1[1])
	//  grad_v0(b) = (v2[2] - v1[2], 0, v1[0] - v2[0])
	//  grad_v0(c) = (v1[1] - v2[1], v2[0] - v1[0], 0)
	g0[0] = (                      b * (v2[2] - v1[2]) + c * (v1[1] - v2[1])) * denom;
	g0[1] = (a * (v1[2] - v2[2]) +                       c * (v2[0] - v1[0])) * denom;
	g0[2] = (a * (v2[1] - v1[1]) + b * (v1[0] - v2[0])                      ) * denom;

	//  grad_v1(a) = (0, v2[2] - v0[2], v0[1] - v2[1])
	//  grad_v1(b) = (v0[2] - v2[2], 0, v2[0] - v0[0])
	//  grad_v1(c) = (v2[1] - v0[1], v0[0] - v2[0], 0)
	g1[0] = (                      b * v02v22          + c * (v2[1] - v0[1])) * denom;
	g1[1] = (a * (v2[2] - v0[2]) +                       c * v00v20         ) * denom;
	g1[2] = (a * v01v21          + b * (v2[0] - v0[0])                      ) * denom;

	//  grad_v2(a) = (0, v0[2] - v1[2], v1[1] - v0[1])
	//  grad_v2(b) = (v1[2] - v0[2], 0, v0[0] - v1[0])
	//  grad_v2(c) = (v0[1] - v1[1], v1[0] - v0[0], 0)
	g2[0] = (                      b * (v1[2] - v0[2]) + c * v01v11) * denom;
	g2[1] = (a * (v02v12) +                       c * (v1[0] - v0[0])) * denom;
	g2[2] = (a * (v1[1] - v0[1]) + b * v00v10                      ) * denom;
}

template<class T>
inline T signed_volume(const Array<T, 3> & v0, const Array<T, 3> & v1, const Array<T, 3> & v2)
{
	// Compute the signed volume of the tetrahedron spanned by the origin, v0, v1, and v2
	return dot(v0, crossProduct(v1, v2)) / 6.;
}

template<class T>
inline void grad_signed_volume(const Array<T, 3> & v0, const Array<T, 3> & v1, const Array<T, 3> & v2,
	Array<T, 3> & grad_v0, Array<T, 3> & grad_v1, Array<T, 3> & grad_v2)
{
	// Use that v0.(v1 x v2) = v2.(v0 x v1) = v1 . (v2 x v0)
	grad_v0 = crossProduct(v1, v2) / 6.;
	grad_v1 = crossProduct(v2, v0) / 6.;
	grad_v2 = crossProduct(v0, v1) / 6.;
}

// Sine and cosine of the angle between the normals of the two triangles
// spanned by (v0, v1, v2) and (v0, v2, v3)
template<class T>
inline void cos_sin_of_angle_between_pair(const Array<T, 3> & v0, const Array<T, 3> & v1,
		const Array<T, 3> & v2, const Array<T, 3> & v3, T & cos_theta, T & sin_theta)
{
	// Compute normals
	Array<T, 3> n1, n2;
	normal(v0, v1, v2, n1);
	normal(v0, v2, v3, n2);

	cos_theta = dot(n1, n2);
	sin_theta = ((dot((n1-n2), v1 - v3) > 0) ? (T)1. : (T)-1.) *
			norm(crossProduct(n1, n2));
}

template<class T>
inline void grad_angle_between_pair(const Array<T, 3> & v0, const Array<T, 3> & v1, const Array<T, 3> & v2, const Array<T, 3> & v3,
		T & cos_angle, T & sin_angle,
		Array<T, 3> & g0, Array<T, 3> & g1, Array<T, 3> & g2, Array<T, 3> & g3)
{
	// Edges
	const Array<T, 3> e0 = v1 - v0;
	const Array<T, 3> e1 = v2 - v0;
	const Array<T, 3> e2 = v3 - v0;
	const T e1_norm = norm(e1);
	const Array<T, 3> e1u = e1 / e1_norm;

	// Normals: n1 = e0 x e1, n2 = e1 x e2
	const Array<T, 3> n1 = crossProduct(e0, e1);
	const Array<T, 3> n2 = crossProduct(e1, e2);
	T n1_norm = norm(n1), n2_norm = norm(n2);
	Array<T, 3> n1u = n1 / n1_norm, n2u = n2 / n2_norm;

	// Compute sin and cos
	cos_angle = dot(n1u, n2u);
	sin_angle = ((dot((n1-n2), v1 - v3) > 0) ? (T)1. : (T)-1.) * norm(crossProduct(n1u, n2u));

	g0 = (dot(e1u, e1 - e2) / n2_norm) * n2u + (dot(e1u, e1 - e0) / n1_norm) * n1u;
	g1 = (-e1_norm / n1_norm) * n1u;
	g2 = (dot(e1u, e2) / n2_norm) * n2u + (dot(e1u, e0) / n1_norm) * n1u;
	g3 = (-e1_norm / n2_norm) * n2u;
}

/*
template<class T>
void cos_sin_grad_angle_between_pair(
		const Array<T, 3> & v0, const Array<T, 3> & v1, const Array<T, 3> & v2, const Array<T, 3> & v3,
		T & cos_angle, T & sin_angle,
		Array<Array<T, 3>, 4> & grad_cos, Array<Array<T, 3>, 4> & grad_sin)
{
	T sin_sign = (dot((crossProduct(v1-v0, v2-v0)-crossProduct(v2-v0, v3-v0)), (v0 + v1 + v2) - (v0 + v2 + v3)) > 0) ? 1 : -1;

	// Code generated from Matlab
	const T t2 = v0[0]-v2[0];
	const T t3 = v0[1]-v2[1];
	const T t4 = v0[0]-v1[0];
	const T t6 = v0[1]-v1[1];
	const T t9 = v0[2]-v2[2];
	const T t10 = v0[2]-v1[2];
	const T t12 = v0[1]-v3[1];
	const T t14 = v0[0]-v3[0];
	const T t17 = v0[2]-v3[2];
	const T t68 = v1[0]-v2[0];
	const T t72 = v2[0]-v3[0];
	const T t45 = v1[1]-v2[1];
	const T t46 = v1[2]-v2[2];
	const T t51 = v2[1]-v3[1];
	const T t52 = v2[2]-v3[2];
	const T t7 = t3*t4 - t2*t6;
	const T t13 = t2*t12;
	const T t29 = t3*t14;
	const T t15 = t13-t29;
	const T t18 = t3*t17 - t9*t12;
	const T t25 = t4*t9 - t2*t10;
	const T t28 = t6*t9 - t3*t10;
	const T t31 = t2*t17 - t9*t14;
	const T t37 = t7*t7 + t25*t25 + t28*t28;
	const T t38 = 1.0/std::sqrt(t37);
	const T t40 = t15*t15 + t18*t18 + t31*t31;
	const T t41 = 1.0/std::sqrt(t40);
	const T t3841 = t38*t41;
	const T t42 = t3841 * (t7*t31 - t15*t25);
	const T t43 = t3841 * (t7*t18 - t15*t28);
	const T t44 = t3841 * (t18*t25 - t28*t31);
	const T t4150 = t41/std::pow(t37,3.0/2.0);
	const T t3856 = t38/std::pow(t40,3.0/2.0);
	const T t67 = t18*t3841*t46;
	const T t76 = 1.0/std::sqrt(t42*t42 + t43*t43 + t44*t44);
	const T t77 = t25*t3841*t52;
	const T t78 = t31*t3841*t46;
	const T t79 = t7*t3841*t52;
	const T t80 = t25*t3841*t72;
	const T t81 = t18*t3841*t45;
	const T t88 = t7*t3841*t51;
	const T t89 = t15*t3841*t45;
	const T t90 = t28*t3841*t72;
	const T t91 = t31*t3841*t45;
	const T t95 = t9*t18*t3841;
	const T t97 = (t2*t7 - t9*t28)*2.0;
	const T t99 = t9*t15*t3841;
	const T t100 = t2*t31*t3841;
	const T t101 = t3*t18*t3841;
	const T t105 = t2*t15*t3841;
	const T t106 = t2*t18*t3841;
	const T t113 = t17*t28*t3841;
	const T t114 = t10*t18*t3841;
	const T t55 = (t15*t51 + t31*t52)*2.0;
	const T t70 = (t28*t46 - t7*t68)*2.0;
	const T t74 = (t18*t52 - t15*t72)*2.0;
	const T t415094  = t4150*(t3*t7 + t9*t25)*2.0;
	const T t415049  = t4150 * (t7*t45 + t25*t46)*2.0;
	const T t415084  = t4150 * (t25*t68 + t28*t45)*2.0;
	const T t385687  = t3856 * (t31*t72 + t18*t51)*2.0;
	const T t4150104 = t4150 * (t2*t25 + t3*t28)*2.0;
	const T t4150109 = t4150 * (t6*t7 + t10*t25)*2.0;
	const T t3856112 = t3856 * (t12*t15 + t17*t31)*2.0;
	const T t4150116 = t4150 * (t4*t7 - t10*t28)*2.0;
	const T t3856119 = t3856 * (t14*t15 - t17*t18)*2.0;
	const T t4150130 = t4150 * (t4*t25 + t6*t28) * 2.0;
	const T t3856133 = t3856 * (t14*t31 + t12*t18) * 2.0;
	const T t3856140 = t3856 * (t3*t15 + t9*t31) * 2.0;
	const T t3856143 = t3856 * (t2*t15 - t9*t18) * 2.0;
	const T t3856151 = t3856 * (t2*t31 + t3*t18) * 2.0;
	const T t122 = t7*t17*t3841;
	const T t123 = t10*t15*t3841;
	const T t124 = t14*t25*t3841;
	const T t125 = t12*t28*t3841;
	const T t126 = t4*t31*t3841;
	const T t127 = t6*t18*t3841;
	const T t134 = t7*t14*t3841;
	const T t135 = t6*t15*t3841;
	const T t136 = t14*t28*t3841;
	const T t137 = t4*t18*t3841;
	const T t141 = t9*t28*t3841;
	const T t145 = t9*t25*t3841;
	const T t146 = t7*t9*t3841;
	const T t147 = t2*t25*t3841;
	const T t148 = t3*t28*t3841;
	const T t152 = t3*t7*t3841;
	const T t153 = t2*t28*t3841;

	cos_angle = (t18*t28 + t25*t31 + t7*t15)*t3841;
	sin_angle = sin_sign/t76;

	grad_cos[0][0] = t77 + t78 + t88 + t89 - t415049*0.5*(t7*t15 + t18*t28 + t25*t31) - t3856*t55*0.5*(t7*t15 + t18*t28 + t25*t31);
	grad_sin[0][0] = sin_sign*t76*(t42*(t79+t91-t15*t3841*t46-t25*t3841*t51-t7*t31*t415049*0.5+t15*t25*t415049*0.5-t7*t31*t3856*t55*0.5+t15*t25*t3856*t55*0.5)*2.0+t43*(t81-t28*t3841*t51-t7*t18*t415049*0.5-t7*t18*t3856*t55*0.5+t15*t28*t415049*0.5+t15*t28*t3856*t55*0.5)*2.0+t44*(t67-t28*t3841*t52-t18*t25*t415049*0.5-t18*t25*t3856*t55*0.5+t28*t31*t415049*0.5+t28*t31*t3856*t55*0.5)*2.0)*0.5;
	grad_cos[0][1] = t67-t7*t3841*t72+t28*t3841*t52-t15*t3841*t68-t7*t15*t4150*t70*0.5-t7*t15*t3856*t74*0.5-t18*t28*t4150*t70*0.5-t18*t28*t3856*t74*0.5-t25*t31*t4150*t70*0.5-t25*t31*t3856*t74*0.5;
	grad_sin[0][1] = sin_sign*t76*(t43*(t79+t90-t15*t3841*t46-t18*t3841*t68-t7*t18*t4150*t70*0.5-t7*t18*t3856*t74*0.5+t15*t28*t4150*t70*0.5+t15*t28*t3856*t74*0.5)*2.0+t44*(t77-t78-t18*t25*t4150*t70*0.5-t18*t25*t3856*t74*0.5+t28*t31*t4150*t70*0.5+t28*t31*t3856*t74*0.5)*2.0+t42*(t80-t31*t3841*t68-t7*t31*t4150*t70*0.5-t7*t31*t3856*t74*0.5+t15*t25*t4150*t70*0.5+t15*t25*t3856*t74*0.5)*2.0)*0.5;
	grad_cos[0][2] = -t80-t81-t28*t3841*t51-t31*t3841*t68+t7*t15*t415084*0.5+t7*t15*t385687*0.5+t18*t28*t415084*0.5+t18*t28*t385687*0.5+t25*t31*t415084*0.5+t25*t31*t385687*0.5;
	grad_sin[0][2] = sin_sign*t76*(t44*(t90+t91-t25*t3841*t51-t18*t3841*t68+t18*t25*t415084*0.5+t18*t25*t385687*0.5-t28*t31*t415084*0.5-t28*t31*t385687*0.5)*-2.0+t43*(t88-t89-t7*t18*t415084*0.5-t7*t18*t385687*0.5+t15*t28*t415084*0.5+t15*t28*t385687*0.5)*2.0+t42*(t7*t3841*t72-t15*t3841*t68-t7*t31*t415084*0.5+t15*t25*t415084*0.5-t7*t31*t385687*0.5+t15*t25*t385687*0.5)*2.0)*(-1.0/2.0);
	grad_cos[1][0] = -t3*t15*t3841 - t9*t31*t3841 + t415094*0.5*(t7*t15 + t18*t28 + t25*t31);
	grad_sin[1][0] = sin_sign*t76*(t42*(t99-t3*t31*t3841+t7*t31*t415094*0.5-t15*t25*t415094*0.5)*-2.0+t43*(t101-t7*t18*t415094*0.5+t15*t28*t415094*0.5)*2.0+t44*(t95-t18*t25*t415094*0.5+t28*t31*t415094*0.5)*2.0)*(-0.5);
	grad_cos[1][1] = -t95 + t105 - t4150*t97*0.5*(t7*t15 + t18*t28 + t25*t31);
	grad_sin[1][1] = sin_sign*t76*(t43*(t99+t106-t7*t18*t4150*t97*0.5+t15*t28*t4150*t97*0.5)*2.0+t44*(t9*t31*t3841-t18*t25*t4150*t97*0.5+t28*t31*t4150*t97*0.5)*2.0+t42*(t100-t7*t31*t4150*t97*0.5+t15*t25*t4150*t97*0.5)*2.0)*0.5;
	grad_cos[1][2] = t100 + t101 - t4150104*0.5*(t7*t15 + t18*t28 + t25*t31);
	grad_sin[1][2] = sin_sign*t76*(t44*(t106-t3*t31*t3841 - t18*t25*t4150104*0.5 + t28*t31*t4150104*0.5)*-2.0 + t43*(t3*t15*t3841 + t7*t18*t4150104*0.5 - t15*t28*t4150104*0.5)*2.0 + t42*(t105 + t7*t31*t4150104*0.5 - t15*t25*t4150104*0.5)*2.0)*(-0.5);
	grad_cos[2][0] = t3841*(t6*t15 + t10*t31 - t7*t12 - t17*t25) - t4150109*0.5*(t7*t15 + t18*t28 + t25*t31) + t3856112*0.5*(t7*t15 + t18*t28 + t25*t31);
	grad_sin[2][0] = sin_sign*t76*(t42*(t122+t123-t6*t31*t3841-t12*t25*t3841+t7*t31*t4150109*0.5-t15*t25*t4150109*0.5-t7*t31*t3856112*0.5+t15*t25*t3856112*0.5)*-2.0+t43*(t125+t127-t7*t18*t4150109*0.5+t7*t18*t3856112*0.5+t15*t28*t4150109*0.5-t15*t28*t3856112*0.5)*2.0+t44*(t113+t114-t18*t25*t4150109*0.5+t18*t25*t3856112*0.5+t28*t31*t4150109*0.5-t28*t31*t3856112*0.5)*2.0)*0.5;
	grad_cos[2][1] = -t113+t114+t134-t4*t15*t3841+t7*t15*t4150116*0.5-t7*t15*t3856119*0.5+t18*t28*t4150116*0.5-t18*t28*t3856119*0.5+t25*t31*t4150116*0.5-t25*t31*t3856119*0.5;
	grad_sin[2][1] = sin_sign*t76*(t44*(t10*t31*t3841+t17*t25*t3841-t18*t25*t4150116*0.5+t18*t25*t3856119*0.5+t28*t31*t4150116*0.5-t28*t31*t3856119*0.5)*2.0+t43*(t122+t123+t136+t137-t7*t18*t4150116*0.5+t7*t18*t3856119*0.5+t15*t28*t4150116*0.5-t15*t28*t3856119*0.5)*2.0+t42*(t124+t126-t7*t31*t4150116*0.5+t15*t25*t4150116*0.5+t7*t31*t3856119*0.5-t15*t25*t3856119*0.5)*2.0)*(-1.0/2.0);
	grad_cos[2][2] = t124+t125-t126-t127+t7*t15*t4150130*0.5-t7*t15*t3856133*0.5+t18*t28*t4150130*0.5-t18*t28*t3856133*0.5+t25*t31*t4150130*0.5-t25*t31*t3856133*0.5;
	grad_sin[2][2] = sin_sign*t76*(t44*(t136+t137-t6*t31*t3841-t12*t25*t3841-t18*t25*t4150130*0.5+t18*t25*t3856133*0.5+t28*t31*t4150130*0.5-t28*t31*t3856133*0.5)*-2.0+t43*(t135+t7*t12*t3841+t7*t18*t4150130*0.5-t7*t18*t3856133*0.5-t15*t28*t4150130*0.5+t15*t28*t3856133*0.5)*2.0+t42*(t134+t4*t15*t3841+t7*t31*t4150130*0.5-t15*t25*t4150130*0.5-t7*t31*t3856133*0.5+t15*t25*t3856133*0.5)*2.0)*0.5;
	grad_cos[3][0] = t145 + t152 - t3856140*0.5*(t7*t15 + t18*t28 + t25*t31);
	grad_sin[3][0] = sin_sign*t76*(t42*(t146-t3*t25*t3841-t7*t31*t3856140*0.5+t15*t25*t3856140*0.5)*-2.0+t43*(t148+t7*t18*t3856140*0.5-t15*t28*t3856140*0.5)*2.0+t44*(t141+t18*t25*t3856140*0.5-t28*t31*t3856140*0.5)*2.0)*(-1.0/2.0);
	grad_cos[3][1] = t141 - t2*t7*t3841 + t3856143*0.5*(t7*t15 + t18*t28 + t25*t31);
	grad_sin[3][1] = sin_sign*t76*(t43*(t146+t153+t7*t18*t3856143*0.5-t15*t28*t3856143*0.5)*2.0+t42*(t147+t7*t31*t3856143*0.5-t15*t25*t3856143*0.5)*2.0+t44*(t145+t18*t25*t3856143*0.5-t28*t31*t3856143*0.5)*2.0)*0.5;
	grad_cos[3][2] = -t147 - t148 + t3856151*0.5*(t7*t15 + t18*t28 + t25*t31);
	grad_sin[3][2] = sin_sign*t76*(t44*(t153-t3*t25*t3841+t18*t25*t3856151*0.5-t28*t31*t3856151*0.5)*-2.0+t42*(t2*t7*t3841-t7*t31*t3856151*0.5+t15*t25*t3856151*0.5)*2.0+t43*(t152-t7*t18*t3856151*0.5+t15*t28*t3856151*0.5)*2.0)*(-1.0/2.0);
}*/

/*
 * Green strain tensor for a triangular element
 */
template<class T>
void green_strain_tensor(const Array<T, 3> & v0, const Array<T, 3> & v1, const Array<T, 3> & v2, Array<Array<T, 2>, 2> & D)
{

}

template<class T>
void grad_green_strain_tensor(const Array<T, 3> & v0, const Array<T, 3> & v1, const Array<T, 3> & v2,
		Array<Array<T, 3>, 3> & D11, Array<Array<T, 3>, 3> & D12, Array<Array<T, 3>, 3> & D22)
{

}

template<class T>
void bending_strain_tensor(const Array<Array<T, 3>, 6> & v, const Array<Array<T, 3>, 6> & v0, Array<Array<T, 2>, 2> & ret)
{
	ret[0][0] = norm(v[3] - 2.*v[2] + v[4]) - norm(v0[3] - 2.*v0[2] + v0[4]);
	ret[0][1] = norm(v[0] - v[1] - v[2] + v[3]) - norm(v0[0] - v0[1] - v0[2] + v0[3]);
	ret[1][0] = ret[0][1];
	ret[1][1] = norm(v[3] - 2.*v[1] + v[5]) - norm(v0[3] - 2.*v0[1] + v0[5]);
}

template<class T>
void grad_bending_strain_tensor(const Array<Array<T, 3>, 6> & v, Array<Array<T, 3>, 6> & K11, Array<Array<T, 3>, 6> & K12, Array<Array<T, 3>, 6> & K22)
{
	// Gradient of K11
	Array<T, 3> tmp = v[3] - 2*v[2] + v[4];
	T tmp_norm = norm(tmp);
    K11[0].resetToZero();
    K11[1].resetToZero();
    K11[2] = -2.*tmp / tmp_norm;
    K11[3] = tmp / tmp_norm;
	K11[4] = tmp / tmp_norm;
	K11[5].resetToZero();

	// Gradient of K12
	tmp = v[0] - v[1] - v[2] + v[3];
	tmp_norm = norm(tmp);
	K12[0] = tmp / tmp_norm;
	K12[1] = -tmp / tmp_norm;
	K12[2] = -tmp / tmp_norm;
	K12[3] = tmp / tmp_norm;
	K12[4].resetToZero();
	K12[5].resetToZero();

	// Gradient of K22
	tmp = v[3] - 2.*v[1] + v[5];
	tmp_norm = norm(tmp);
	K22[0].resetToZero();
	K22[1] = -2.*tmp / tmp_norm;
	K22[2].resetToZero();
	K22[3] = tmp / tmp_norm;
	K22[4].resetToZero();
	K22[5] = tmp / tmp_norm;
}

} /* namespace tri */

} /* namespace fsi */

} /* namespace plb */



#endif /* TRIANGLEUTILS_H_ */
