/*
 * Grid.h
 *
 *  Created on: Mar 26, 2014
 *      Author: niber
 */

#ifndef GRID_H_
#define GRID_H_
#include "geometry.h"

namespace plb {

namespace fsi {

template<class T, class U>
struct Cell {
	std::vector<U> nodes;
	geo::Rect<T> domain;
};

template<class T, class U>
class Grid {
public:
	typedef Cell<T, U> CellType;

	Grid();
	Grid(const geo::Rect<T> domain_, T min_cell_size_);
	CellType & get_cell(plint i, plint j, plint k);
	void insert(const U & node, const Array<T, 3> & pos);
	Dot3D get_index(const Array<T, 3> &) const;
	void clear();
	void repartition();

	void set_min_cell_size(T new_size) { min_cell_size = new_size; }
	T get_min_cell_size() const { return min_cell_size; }
	void set_domain(const geo::Rect<T> & new_domain);
	geo::Rect<T> get_domain() const { return domain; }

	Box3D get_bulk_indices() const {
		return Box3D(
				1, num_cells.x-2,
				1, num_cells.y-2,
				1, num_cells.z-2);
	}

	pluint cell_count_x() const { return num_cells.x; }
	pluint cell_count_y() const { return num_cells.y; }
	pluint cell_count_z() const { return num_cells.z; }
	pluint cell_count() const { return cells.size(); }

private:
	void copy_nodes_from_cell_to_cell(const CellType &, CellType &);

	T x0, y0, z0;
	T min_cell_size;
	std::vector<CellType> cells;
	Array<T, 3> cell_size;
	Dot3D num_cells;
	geo::Rect<T> domain;
};

}

}

#endif /* GRID_H_ */
