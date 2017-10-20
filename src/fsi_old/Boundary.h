#ifndef BOUNDARY_H_
#define BOUNDARY_H_
#include "geometry.h"

namespace plb {
namespace fsi {

template<class T>
class Boundary {
public:
	virtual ~Boundary() { }
	virtual bool contains(const Array<T, 3> &, T margin = 0) const = 0;
	virtual bool does_intersect(const geo::Rect<T> &, T margin = 0) const = 0;
	virtual bool distance_to_boundary_less_than(const Array<T, 3> &, T) const = 0;
	virtual T distance_to_boundary(const Array<T, 3> &) const = 0;
	virtual T distance_to_boundary(const Array<T, 3> &, const Array<T, 3> &) const = 0;
	virtual Array<T, 3> get_normal(const Array<T, 3> &) const = 0;
};

/*
 * Pipe boundary aligned with the x-direction
 */
template<class T>
class PipeBoundary : public Boundary<T> {
public:
	PipeBoundary(T y0, T z0, T radius)
	: y0(y0), z0(z0), radius(radius)
	{ }

	virtual ~PipeBoundary() { }

	virtual bool contains(const Array<T,3> & pos, T margin = 0) const
	{
		return (util::sqr(pos[1]-y0) + util::sqr(pos[2]-z0)) <= util::sqr(radius - margin);
	}

	virtual bool distance_to_boundary_less_than(const Array<T, 3> & pos, T dist) const
	{
		return (util::sqr(pos[1]-y0) + util::sqr(pos[2]-z0)) > util::sqr(radius-dist);
	}

	virtual T distance_to_boundary(const Array<T, 3> & pos) const
	{
		return radius - std::sqrt(util::sqr(pos[1]-y0) + util::sqr(pos[2]-z0));
	}

	/**
	 * distance_to_boundary
	 *   Params: pos - position to trace from
	 *           dir - direction in which the distance to the boundary is traced
	 *   Returns the smallest distance (in absolute sense) to the boundary
	 */
	virtual T distance_to_boundary(const Array<T, 3> & pos, const Array<T, 3> & dir) const
	{
		const T dy = pos[1] - y0;
		const T dz = pos[2] - z0;
		const T nnrm_sqr = util::sqr(dir[1]) + util::sqr(dir[2]);

		const T a = (dir[1]*dy+ dir[2]*dz) / nnrm_sqr;
		const T b = a*a - (dy*dy + dz*dz - radius*radius)/nnrm_sqr;

		// No solution
		if(b <= 0.0001)
			return std::numeric_limits<T>::max();

		// Find the smallest solution
		const T c = std::sqrt(b);

		const T dist1 = -a + c;
		const T dist2 = -a - c;

		return std::abs(dist1) < std::abs(dist2) ? dist1 : dist2;
	}

	/*virtual void get_mirror_point(const Array<T, 3> & p, Array<T, 3> & res) const
	{
		const T factor = (T) 2 * (radius / std::sqrt(util::sqr(p[1]-y0) + util::sqr(p[2]-z0)) - (T)1);
		res[0] = p[0];
		res[1] = p[1] + factor * (p[1] - y0);
		res[2] = p[2] + factor * (p[2] - z0);
	}*/

	virtual bool does_intersect(const geo::Rect<T> & r, T width = 0) const
	{
		// Calculate the distance between the cylinder center and the closest point in the rectangle
		T tmp = (util::sqr(y0 - geo::clamp(y0, r.y0, r.y1)) +
				util::sqr(z0 - geo::clamp(z0, r.z0, r.z1)));
		// Test if the rectangle intersects the cylinder
		if(tmp <= util::sqr(radius)) {
			// Check if the rectangle is fully contained in the cylinder with radius = radius - width
			// In that case, the rectangle does not intersect the boundary
			return ! (contains(Array<T, 3>((T)0, r.y0, r.z0), width) &&
					  contains(Array<T, 3>((T)0, r.y0, r.z1), width) &&
					  contains(Array<T, 3>((T)0, r.y1, r.z0), width) &&
					  contains(Array<T, 3>((T)0, r.y1, r.z1), width));
		} else {
			return false;
		}
	}

	virtual Array<T, 3> get_normal(const Array<T, 3> & pos) const
	{
		Array<T, 3> d(0, pos[1]-y0, pos[2]-z0);
		return d / norm(d);
	}

	T y0, z0, radius;
};

/*
 * Couette flow with walls in the xz-plane
 */
/*template<class T>
class CouetteBoundary : public Boundary<T> {
public:
	CouetteBoundary(T y0_, T y1_) : y0(y0_), y1(y1_) { }

	virtual ~CouetteBoundary() { }

	virtual bool contains(const Array<T, 3> & pos, T margin = 0) const
	{
		return pos[1] > (y0+margin) && pos[1] < (y1-margin);
	}

	virtual bool does_intersect(const geo::Rect<T> &, T margin = 0) const
	{

	}

	virtual bool distance_to_boundary_less_than(const Array<T, 3> &, T) const
	{

	}

	virtual void get_mirror_point(const Array<T, 3> & pos, Array<T, 3> & ret) const
	{
		ret[0] = pos[0];
		ret[2] = pos[2];
		if(std::abs(pos[1] - y0) < std::abs(pos[1] - y1)) {
			ret[1] =
		} else {

		}
	}
private:
	T y0, y1;
};*/

}

}




#endif /* BOUNDARY_H_ */
