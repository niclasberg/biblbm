/*
 * TriangleHash.hh
 *
 *  Created on: Jun 22, 2015
 *      Author: niber
 */

#ifndef TRIANGLEHASH_HH_
#define TRIANGLEHASH_HH_
#include "TriangleHash.h"
#include "Periodicity.h"

namespace plb {

namespace fsi {

template<class T, class TriangleType, class VertexType>
TriangleHash<T, TriangleType, VertexType>::TriangleHash(const std::vector<TriangleType> & triangles, const std::vector<VertexType> & vertices)
: triangles_(triangles), vertices_(vertices), offset_(0, 0, 0), bounding_box_(0,0,0,0,0,0), need_rehash_(true)
{

}

template<class T, class TriangleType, class VertexType>
void TriangleHash<T, TriangleType, VertexType>::rehash()
{
	// Recompute bounding box
	bounding_box_.x0 = std::numeric_limits<plint>::max();
	bounding_box_.y0 = std::numeric_limits<plint>::max();
	bounding_box_.z0 = std::numeric_limits<plint>::max();
	bounding_box_.x1 = -std::numeric_limits<plint>::max();
	bounding_box_.y1 = -std::numeric_limits<plint>::max();
	bounding_box_.z1 = -std::numeric_limits<plint>::max();

	for(typename std::vector<VertexType>::const_iterator it = vertices_.begin(); it != vertices_.end(); ++it) {
		bounding_box_.x0 = std::min(bounding_box_.x0, (plint)std::floor(it->pos[0]));
		bounding_box_.y0 = std::min(bounding_box_.y0, (plint)std::floor(it->pos[1]));
		bounding_box_.z0 = std::min(bounding_box_.z0, (plint)std::floor(it->pos[2]));
		bounding_box_.x1 = std::max(bounding_box_.x1, (plint)std::ceil(it->pos[0]));
		bounding_box_.y1 = std::max(bounding_box_.y1, (plint)std::ceil(it->pos[1]));
		bounding_box_.z1 = std::max(bounding_box_.z1, (plint)std::ceil(it->pos[2]));
	}

	offset_.x = bounding_box_.x0;
	offset_.y = bounding_box_.y0;
	offset_.z = bounding_box_.z0;

	bounding_box_ = bounding_box_.shift(-offset_.x, -offset_.y, -offset_.z);

	// Resize the hash container so that it can hold all cells
	triangle_hash_.resize(num_cells_x()*num_cells_y()*num_cells_z());

	clear();

	// Hash the triangles
	for(plint i_tri = 0; i_tri < triangles_.size(); ++i_tri) {
		const TriangleType & tri = triangles_[i_tri];
		Box3D bb = compute_bounding_box(tri);

		for(plint i = bb.x0; i <= bb.x1; ++i)
			for(plint j = bb.y0; j <= bb.y1; ++j)
				for(plint k = bb.z0; k <= bb.z1; ++k) {
					this->get(i, j, k).triangles.push_back(&tri);
				}
	}
}

template<class T, class TriangleType, class VertexType>
void TriangleHash<T, TriangleType, VertexType>::clear()
{
	for(plint i = 0; i < triangle_hash_.size(); ++i)
		triangle_hash_[i].triangles.clear();
}

template<class T, class TriangleType, class VertexType>
Box3D TriangleHash<T, TriangleType, VertexType>::compute_bounding_box(const TriangleType & tri)
{
	PLB_PRECONDITION(tri.i0 < vertices_.size() && tri.i1 < vertices_.size() && tri.i2 < vertices_.size())
	return Box3D(
		std::floor(std::min(vertices_[tri.i0].pos[0], std::min(vertices_[tri.i1].pos[0], vertices_[tri.i2].pos[0]))) - offset_.x,
		std::ceil(std::max(vertices_[tri.i0].pos[0], std::max(vertices_[tri.i1].pos[0], vertices_[tri.i2].pos[0]))) - offset_.x,
		std::floor(std::min(vertices_[tri.i0].pos[1], std::min(vertices_[tri.i1].pos[1], vertices_[tri.i2].pos[1]))) - offset_.y,
		std::ceil(std::max(vertices_[tri.i0].pos[1], std::max(vertices_[tri.i1].pos[1], vertices_[tri.i2].pos[1]))) - offset_.y,
		std::floor(std::min(vertices_[tri.i0].pos[2], std::min(vertices_[tri.i1].pos[2], vertices_[tri.i2].pos[2]))) - offset_.z,
		std::ceil(std::max(vertices_[tri.i0].pos[2], std::max(vertices_[tri.i1].pos[2], vertices_[tri.i2].pos[2]))) - offset_.z
	);
}

template<class T, class TriangleType, class VertexType>
void TriangleHash<T, TriangleType, VertexType>::voxelize()
{
	// Do ray-tracing in the z-direction, starting from outside the bounding domain.
	// if the ray has crossed zero or an even number of triangles, the nodes crossed by the
	// ray are outside of the object, and vice versa.
	T dist;

	// Use std::sets for the candidates. This ensure that all elements will be unique
	// and duplicate ray-triangle intersection tests are avoided.
	std::set<const TriangleType *> candidates;
	std::set<const TriangleType *> candidates_next;

	for(plint i = bounding_box_.x0; i <= bounding_box_.x1; ++i) {
		for(plint j = bounding_box_.y0; j <= bounding_box_.y1; ++j) {
			bool is_inside = false;
			candidates.clear();
			candidates_next.clear();

			for(plint k = bounding_box_.z0; k <= bounding_box_.z1; ++k) {
				//std::cout << get_ind(i, j, k) << " / " << triangle_hash_.size() << std::endl;
				List & list = get(i, j, k);

				// Find unique candidates for overlap
				candidates = candidates_next;
				for(plint it = 0; it < list.triangles.size(); ++it)
					candidates.insert(list.triangles[it]);
				//std::cout << candidates.size() << std::endl;

				candidates_next.clear();

				for(typename std::set<const TriangleType *>::const_iterator it = candidates.begin(); it != candidates.end(); ++it) {
					if(tri_ray_intersect(*(*it), Array<T, 3>(offset_.x+i, offset_.y+j, offset_.z+k-1), Array<T, 3>(0, 0, 1), dist)) {
						//std::cout << dist << std::endl;
						if(dist >= 0) {
							if(dist < 1)
								is_inside = !is_inside;
							else
								// Overlap after the cell, store for next iteration
								candidates_next.insert(*it);
						}
					}
				}

				// Insert value
				list.is_inside = is_inside;
			}
		}
	}
}

template<class T, class TriangleType, class VertexType>
bool TriangleHash<T, TriangleType, VertexType>::tri_ray_intersect(const TriangleType & tri, const Array<T, 3> & pos, const Array<T, 3> & dir, T & dist) const
{
	static T epsilon = 0.00001;

	//std::cout << vertices_[tri.i0].pos[2] << ", " << pos[2] << std::endl;

	const Array<T, 3> e1 = vertices_[tri.i1].pos-vertices_[tri.i0].pos;
	const Array<T, 3> e2 = vertices_[tri.i2].pos-vertices_[tri.i0].pos;
	const Array<T, 3> q  = crossProduct(dir,e2);
	const T a  = dot(e1,q);

	// Check if the vector is parallel to the plane (the intersection is at infinity)
	if (a>-epsilon && a<epsilon)
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

template<class T, class TriangleType, class VertexType>
void TriangleHash<T, TriangleType, VertexType>::to_scalar_field(const Box3D & domain, ScalarField3D<T> & field) const
{
	to_scalar_field(domain, field, NormalArithmetic<T>());
}

template<class T, class TriangleType, class VertexType>
	template<class Arithmetic>
void TriangleHash<T, TriangleType, VertexType>::to_scalar_field(const Box3D & domain, ScalarField3D<T> & field, const Arithmetic & arithmetic) const
{
	// TODO fix the indices here!
	Dot3D field_bb_offset = offset_ - field.getLocation();
	Box3D intersection;
	if(intersect(bounding_box_, domain.shift(-offset_.x, -offset_.y, -offset_.z), intersection)) {
		for(plint i = intersection.x0; i <= intersection.x1; ++i)
			for(plint j = intersection.y0; j <= intersection.y1; ++j)
				for(plint k = intersection.z0; k <= intersection.z1; ++k) {
					if(get(i, j, k).is_inside)
						field.get(i+field_bb_offset.x, j+field_bb_offset.y, k+field_bb_offset.z) = 1;
				}
	}
}

template<class T, class TriangleType, class VertexType>
void TriangleHash<T, TriangleType, VertexType>::to_dot_list(const Box3D & domain, std::vector<Dot3D> & inner_pts) const
{
	to_dot_list(domain, inner_pts, NormalArithmetic<T>());
}

template<class T, class TriangleType, class VertexType>
	template<class Arithmetic>
void TriangleHash<T, TriangleType, VertexType>::to_dot_list(const Box3D & domain, std::vector<Dot3D> & inner_pts, const Arithmetic & arithmetic) const
{
	Box3D intersection;
	plint pXmin, pXmax, pYmin, pYmax, pZmin, pZmax;
	pXmin = pXmax = pYmin = pYmax = pZmin = pZmax = 0;
	if(arithmetic.periodic_x()) {
		pXmin = -1;
		pXmax = 1;
	}
	if(arithmetic.periodic_y()) {
		pYmin = -1;
		pYmax = 1;
	}
	if(arithmetic.periodic_z()) {
		pZmin = -1;
		pZmax = 1;
	}

	for(plint sX = pXmin; sX <= pXmax; ++sX)
		for(plint sY = pYmin; sY <= pYmax; ++sY)
			for(plint sZ = pZmin; sZ <= pZmax; ++sZ) {
				// Periodically shifted domain
				Box3D domain2 = domain.shift(sX*arithmetic.get_nx(),
											 sY*arithmetic.get_ny(),
											 sZ*arithmetic.get_nz());

				if(intersect(bounding_box_, domain2.shift(-offset_.x, -offset_.y, -offset_.z), intersection)) {
					for(plint i = intersection.x0; i <= intersection.x1; ++i)
						for(plint j = intersection.y0; j <= intersection.y1; ++j)
							for(plint k = intersection.z0; k <= intersection.z1; ++k) {
								if(get(i, j, k).is_inside)
									inner_pts.push_back(
											Dot3D(i+offset_.x - sX*arithmetic.get_nx(),
												  j+offset_.y - sY*arithmetic.get_ny(),
												  k+offset_.z - sZ*arithmetic.get_nz()));
							}
				}
			}
}

}

}



#endif /* TRIANGLEHASH_HH_ */
