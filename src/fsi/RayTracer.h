/*
 * RayTracer.h
 *
 *  Created on: Jul 23, 2015
 *      Author: niber
 */

#ifndef RAYTRACER_H_
#define RAYTRACER_H_
#include <vector>

namespace plb {

namespace fsi {

template<class T, class TriangleType, class VertexType>
class RayTracer {
public:
	RayTracer(const std::vector<TriangleType> & triangles, const std::vector<VertexType> & vertices)
	: vertices_(vertices), triangles_(triangles)
	{ }

	template<class Container> void find_intersections_xyz(Container &, Container &, Container &) const;
	template<class Container> void find_intersections_z(Container &) const;

private:
	bool tri_ray_intersect(const TriangleType & tri, const Array<T, 3> & pos, const Array<T, 3> & dir, T & dist) const;
	void compute_bounding_box(const TriangleType & tri, Box3D &) const;

private:
	const std::vector<VertexType> & vertices_;
	const std::vector<TriangleType> & triangles_;

};

}

}



#endif /* RAYTRACER_H_ */
