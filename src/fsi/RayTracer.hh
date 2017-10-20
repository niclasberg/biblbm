/*
 * RayTracer.hh
 *
 *  Created on: Jul 23, 2015
 *      Author: niber
 */

#ifndef RAYTRACER_HH_
#define RAYTRACER_HH_
#include "RayTracer.h"
#include "ZBuffer.h"

namespace plb {

namespace fsi {

template<class T, class TriangleType, class VertexType>
	template<class Container>
void RayTracer<T, TriangleType, VertexType>::find_intersections_z(Container & z_buffer) const
{
	T dist;
	Array<T, 3> normal;

	for(plint iTri = 0; iTri < triangles_.size(); ++iTri) {
		const TriangleType & tri = triangles_[iTri];

		const plint iXmin = std::ceil(std::min(vertices_[tri.i0].pos[0], std::min(vertices_[tri.i1].pos[0], vertices_[tri.i2].pos[0])));
		const plint iXmax = std::floor(std::max(vertices_[tri.i0].pos[0], std::max(vertices_[tri.i1].pos[0], vertices_[tri.i2].pos[0])));
		const plint iYmin = std::ceil(std::min(vertices_[tri.i0].pos[1], std::min(vertices_[tri.i1].pos[1], vertices_[tri.i2].pos[1])));
		const plint iYmax = std::floor(std::max(vertices_[tri.i0].pos[1], std::max(vertices_[tri.i1].pos[1], vertices_[tri.i2].pos[1])));

		for(plint i = iXmin; i <= iXmax; ++i)
			for(plint j = iYmin; j <= iYmax; ++j) {
				if(tri_ray_intersect(triangles_[iTri], Array<T, 3>(i, j, 0.), Array<T, 3>(0, 0, 1), dist)) {
					z_buffer.put(i, j, dist);
				}
			}
	}
}

template<class T, class TriangleType, class VertexType>
	template<class Container>
void RayTracer<T, TriangleType, VertexType>::find_intersections_xyz(Container & x_buffer, Container & y_buffer, Container & z_buffer) const
{
	T dist;
	Array<T, 3> normal;
	Array<T, 3> startPoint;
	Array<T, 3> traceDirections[3] = {Array<T, 3>(1., 0, 0), Array<T, 3>(0, 1., 0), Array<T, 3>(0, 0, 1.)};
	Container * buffers[3] = {&x_buffer, &y_buffer, &z_buffer};

	// Ray tracing planes
	plint dirs[3][2] = {{1, 2}, {0, 2}, {0, 1}};

	for(plint iTri = 0; iTri < triangles_.size(); ++iTri) {
		const TriangleType & tri = triangles_[iTri];

		for(plint iDir = 0; iDir < 3; ++iDir) {
			const plint d0 = dirs[iDir][0];
			const plint d1 = dirs[iDir][1];

			const plint i0min = std::ceil(std::min(vertices_[tri.i0].pos[d0], std::min(vertices_[tri.i1].pos[d0], vertices_[tri.i2].pos[d0])));
			const plint i0max = std::floor(std::max(vertices_[tri.i0].pos[d0], std::max(vertices_[tri.i1].pos[d0], vertices_[tri.i2].pos[d0])));
			const plint i1min = std::ceil(std::min(vertices_[tri.i0].pos[d1], std::min(vertices_[tri.i1].pos[d1], vertices_[tri.i2].pos[d1])));
			const plint i1max = std::floor(std::max(vertices_[tri.i0].pos[d1], std::max(vertices_[tri.i1].pos[d1], vertices_[tri.i2].pos[d1])));

			startPoint.resetToZero();
			for(plint i = i0min; i <= i0max; ++i)
				for(plint j = i1min; j <= i1max; ++j) {
					startPoint[d0] = i;
					startPoint[d1] = j;
					if(tri_ray_intersect(triangles_[iTri], startPoint, traceDirections[iDir], dist)) {
						buffers[iDir]->put(i, j, dist);
					}
				}
		}
	}
}

template<class T, class TriangleType, class VertexType>
bool RayTracer<T, TriangleType, VertexType>::tri_ray_intersect(const TriangleType & tri, const Array<T, 3> & pos, const Array<T, 3> & dir, T & dist) const
{
	const Array<T, 3> e1 = vertices_[tri.i1].pos-vertices_[tri.i0].pos;
	const Array<T, 3> e2 = vertices_[tri.i2].pos-vertices_[tri.i0].pos;
	const Array<T, 3> q  = crossProduct(dir,e2);
	const T a  = dot(e1,q);

	// Check if the vector is parallel to the plane (the intersection is at infinity)
	if (std::abs(a) < std::numeric_limits<T>::epsilon())
		return false;

	const T f = 1/a;
	const Array<T, 3> s = pos - vertices_[tri.i0].pos;
	const T u = f*dot(s,q);

	// Check if the intersection is outside of the triangle
	if (u<0.0)
		return false;

	const Array<T, 3> r = crossProduct(s,e1);
	const T v = f*dot(dir,r);

	// Check if the intersection is outside of the triangle
	if(v<0.0 || u+v>1.0)
		return false;

	dist = f*dot(e2,r);
	return true;
}

}

}


#endif /* RAYTRACER_HH_ */
