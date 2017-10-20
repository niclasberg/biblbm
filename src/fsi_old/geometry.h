/*
 * geometry.h
 *
 *  Created on: Mar 4, 2014
 *      Author: niber
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_
#include "core/geometry3D.h"
#include "core/util.h"

namespace plb {

namespace fsi {

namespace geo {

template<class T>
T clamp(const T & val, const T & lower, const T & upper)
{
	if(val >= upper)
		return upper;
	if(val <= lower)
		return lower;
	return val;
}

template<class T>
struct Sphere {
	Sphere() : radius(0), pos(0, 0, 0) { }
	Sphere(T r, const Array<T, 3> & p) : radius(r), pos(p) { }

	bool contains(const Array<T, 3> & p)
	{
		return (util::sqr(p[0] - pos[0]) +
			    util::sqr(p[1] - pos[1]) +
			    util::sqr(p[2] - pos[2]))
				<= util::sqr(radius);
	}

	template<class Arithmetic>
	bool contains(const Array<T, 3> & p, const Arithmetic & a)
	{
		return a.dist_sqr(p, pos) <= util::sqr(radius);
	}

	Sphere & enlarge_inplace(T dr) { radius += dr; return *this;}
	T volume() { return (T) 4 * M_PI * radius * radius * radius / (T) 3; }

	T radius;
	Array<T, 3> pos;
};

template<class T>
struct Rect {
	Rect() : x0(0), x1(0), y0(0), y1(0), z0(0), z1(0) { }
	Rect(T x0_, T x1_, T y0_, T y1_, T z0_, T z1_) : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_) { }
	Rect(const Box3D & b) : x0(b.x0), x1(b.x1), y0(b.y0), y1(b.y1), z0(b.z0), z1(b.z1) { }

	Rect enlarge(T dx, T dy, T dz) {
		Rect ret;
		ret.x0 = x0-dx;
		ret.x1 = x1+dx;
		ret.y0 = y0-dy;
		ret.y1 = y1+dy;
		ret.z0 = z0-dz;
		ret.z1 = z1+dz;
		return ret;
	}

	Rect enlarge(T delta) { return enlarge(delta, delta, delta); }

	Rect & enlarge_inplace(T dx, T dy, T dz)
	{
		x0 -= dx;
		x1 += dx;
		y0 -= dy;
		y1 += dy;
		z0 -= dz;
		z1 += dz;
		return *this;
	}

	void enlarge_inplace(T delta) { enlarge_inplace(delta, delta, delta); }

	bool contains(const Array<T, 3> & p)
	{
		return (p[0] > x0) && (p[0] < x1) &&
			   (p[1] > y0) && (p[1] < y1) &&
			   (p[2] > z0) && (p[2] < z1);
	}

	template<class Arithmetic>
	bool contains(const Array<T, 3> & p, const Arithmetic & a)
	{
		return std::abs(a.dist_x(0.5*(x0+x1), p[0])) < 0.5*(x1-x0) &&
				std::abs(a.dist_y(0.5*(y0+y1), p[1])) < 0.5*(y1-y0) &&
				std::abs(a.dist_z(0.5*(z0+z1), p[2])) < 0.5*(z1-z0);
	}

	bool contains_or_on_boundary(const Array<T, 3> & p)
	{
		return (p[0] >= x0) && (p[0] <= x1) &&
			   (p[1] >= y0) && (p[1] <= y1) &&
			   (p[2] >= z0) && (p[2] <= z1);
	}

	template<class Arithmetic>
	bool contains_or_on_boundary(const Array<T, 3> & p, const Arithmetic & a)
	{
		return std::abs(a.dist_x(0.5*(x0+x1), p[0])) <= 0.5*(x1-x0) &&
				std::abs(a.dist_y(0.5*(y0+y1), p[1])) <= 0.5*(y1-y0) &&
				std::abs(a.dist_z(0.5*(z0+z1), p[2])) <= 0.5*(z1-z0);
	}

	T volume() { return (x1-x0)*(y1-y0)*(z1-z0); }

	T x0, x1, y0, y1, z0, z1;
};

/*
 * Intersection test methods
 */
template<class T>
bool does_intersect(const Rect<T> & r, const Sphere<T> & s)
{
	return (util::sqr(clamp(s.pos[0], r.x0, r.x1) - s.pos[0]) +
			util::sqr(clamp(s.pos[1], r.y0, r.y1) - s.pos[1]) +
			util::sqr(clamp(s.pos[2], r.z0, r.z1) - s.pos[2]))
			<= util::sqr(s.radius);
}

template<class T>
bool does_intersect(const Sphere<T> & s, const Rect<T> & r)
{
	return does_intersect(r, s);
}

template<class T>
bool does_intersect(const Box3D & r, const Sphere<T> & s)
{
	return (util::sqr(clamp(s.pos[0], (T)r.x0, (T)r.x1) - s.pos[0]) +
			util::sqr(clamp(s.pos[1], (T)r.y0, (T)r.y1) - s.pos[1]) +
			util::sqr(clamp(s.pos[2], (T)r.z0, (T)r.z1) - s.pos[2]))
			<= util::sqr(s.radius);
}

template<class T>
bool does_intersect(const Sphere<T> & s, const Box3D & r)
{
	return does_intersect(r, s);
}

template<class T>
bool does_intersect(const Sphere<T> & s, const Sphere<T> & s2)
{
	return (util::sqr(s.pos[0] - s2.pos[0]) +
			util::sqr(s.pos[1] - s2.pos[1]) +
			util::sqr(s.pos[2] - s2.pos[2]))
				<= util::sqr(s.radius + s2.radius);
}

/*template<class T>
bool does_intersect(const CylinderX<T> & c, const Rect<T> & r)
{
	if(c.origin[0] > r.x0 || (c.origin[0] + c.length) > r.x1) {
		return (util::sqr(c.origin[1] - clamp(c.origin[1], r.y0, r.y1)) +
				util::sqr(c.origin[2] - clamp(c.origion[2], r.z0, r.z1)))
				< util::sqr(c.radius);
	}
	return false;
}

template<class T>
bool does_intersect(const Rect<T> & r, const CylinderX<T> & c)
{
	return does_intersect(c, r);
}*/


// Intersection tests with custom arithmetic
template<class T, class Arithmetic>
bool does_intersect(const Rect<T> & r, const Rect<T> & r2, const Arithmetic & a)
{
	// Check if distance between centers is less than the sum of the half-length of the rectangles
	return ! (std::abs(a.dist_x(0.5*(r.x0+r.x1), 0.5*(r2.x0+r2.x1))) > 0.5*(r.x1 - r.x0 + r2.x1 - r2.x0) ||
		      std::abs(a.dist_y(0.5*(r.y0+r.y1), 0.5*(r2.y0+r2.y1))) > 0.5*(r.y1 - r.y0 + r2.y1 - r2.y0) ||
		      std::abs(a.dist_z(0.5*(r.z0+r.z1), 0.5*(r2.z0+r2.z1))) > 0.5*(r.z1 - r.z0 + r2.z1 - r2.z0));
}

template<class T, class Arithmetic>
bool does_intersect(const Rect<T> & r, const Sphere<T> & s, const Arithmetic & a)
{
	// Shift position as close to the rectangle as possible
	Array<T, 3> pos = s.pos;
	a.shift_periodically_to_minimize_distance_to(
			Array<T, 3>(0.5*(r.x0+r.x1), 0.5*(r.y0+r.y1), 0.5*(r.z0+r.z1)),
			pos);

	// Do usual test
	return (util::sqr(clamp(pos[0], r.x0, r.x1) - pos[0]) +
			util::sqr(clamp(pos[1], r.y0, r.y1) - pos[1]) +
			util::sqr(clamp(pos[2], r.z0, r.z1) - pos[2]))
			<= util::sqr(s.radius);
}

template<class T, class Arithmetic>
bool does_intersect(const Sphere<T> & s, const Rect<T> & r, const Arithmetic & a)
{
	return does_intersect(r, s, a);
}

}

}

}



#endif /* GEOMETRY_H_ */
