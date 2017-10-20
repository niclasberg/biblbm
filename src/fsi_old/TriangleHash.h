/*
 * TriangleHash.h
 *
 *  Created on: Jun 18, 2015
 *      Author: niber
 */

#ifndef TRIANGLEHASH_H_
#define TRIANGLEHASH_H_
#include "core/plbDebug.h"

namespace plb {

namespace fsi {

/*
 * VertexType must have a field "pos" with type Array<T, 3>
 * TriangleType must have the fields i0, i1, i2 of integer type
 */
template<class T, class TriangleType, class VertexType>
class TriangleHash {
public:
	struct List {
		std::vector<const TriangleType *> triangles;
		bool is_inside;
	};

	TriangleHash(const std::vector<TriangleType> & triangles, const std::vector<VertexType> & vertices);

	void rehash();
	void voxelize();
	void clear();

	const List & get(plint i, plint j, plint k) const { return triangle_hash_[get_ind(i, j, k)]; }

	// Data extraction
	template<class Arithmetic>
	void to_scalar_field(const Box3D &, ScalarField3D<T> &, const Arithmetic &) const;
	void to_scalar_field(const Box3D &, ScalarField3D<T> &) const;

	template<class Arithmetic>
	void to_dot_list(const Box3D &, std::vector<Dot3D> &, const Arithmetic &) const;
	void to_dot_list(const Box3D &, std::vector<Dot3D> &) const;

private:
	List & get(plint i, plint j, plint k) { return triangle_hash_[get_ind(i, j, k)]; }

	bool tri_ray_intersect(const TriangleType & tri, const Array<T, 3> & pos, const Array<T, 3> & dir, T & dist) const;
	Box3D compute_bounding_box(const TriangleType & tri);
	plint get_ind(plint i, plint j, plint k) const { return k + num_cells_z()*(j + num_cells_y()*i); }
	plint num_cells_x() const { return bounding_box_.x1 - bounding_box_.x0 + 1; }
	plint num_cells_y() const { return bounding_box_.y1 - bounding_box_.y0 + 1; }
	plint num_cells_z() const { return bounding_box_.z1 - bounding_box_.z0 + 1; }

private:

	bool need_rehash_;
	Dot3D offset_;
	Box3D bounding_box_;
	std::vector<List> triangle_hash_;
	const std::vector<VertexType> & vertices_;
	const std::vector<TriangleType> & triangles_;
};

}

}



#endif /* TRIANGLEHASH_H_ */
