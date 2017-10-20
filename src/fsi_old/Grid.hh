/*
 * Grid.hh
 *
 *  Created on: Mar 26, 2014
 *      Author: niber
 */

#ifndef GRID_HH_
#define GRID_HH_

#include "Grid.h"

namespace plb {

namespace fsi {

template<class T, class U>
Grid<T, U>::Grid() : domain(0, 1, 0, 1, 0, 1), min_cell_size(1) { }

template<class T, class U>
Grid<T, U>::Grid(const geo::Rect<T> domain_, T min_cell_size_)
: domain(domain_), min_cell_size(min_cell_size_)
{
	repartition();
}

template<class T, class U>
void Grid<T, U>::repartition()
{
	// Clear current structure
	cells.clear();

	// Determine number of inner cells
	num_cells.x = std::floor((domain.x1-domain.x0) / min_cell_size);
	num_cells.y = std::floor((domain.y1-domain.y0) / min_cell_size);
	num_cells.z = std::floor((domain.z1-domain.z0) / min_cell_size);

	// Set cell sizes
	cell_size[0] = (domain.x1 - domain.x0) / (T) num_cells.x;
	cell_size[1] = (domain.y1 - domain.y0) / (T) num_cells.y;
	cell_size[2] = (domain.z1 - domain.z0) / (T) num_cells.z;

	// Increase the size of the domain to include an envelope
	num_cells.x += 2;
	num_cells.y += 2;
	num_cells.z += 2;

	x0 = domain.x0 - cell_size[0];
	y0 = domain.y0 - cell_size[1];
	z0 = domain.z0 - cell_size[2];

	// Allocate new cells
	cells.resize(num_cells.x * num_cells.y * num_cells.z);

	// Set domains for each cell
	for(plint i = 0; i < cell_count_x(); ++i) {
		T x0 = domain.x0 + (i-1) * cell_size[0];
		for(plint j = 0; j < cell_count_y(); ++j) {
			T y0 = domain.y0 + (j-1) * cell_size[1];
			for(plint k = 0; k < cell_count_z(); ++k) {
				T z0 = domain.z0 + (k-1) * cell_size[2];
				get_cell(i, j, k).domain = geo::Rect<T>(
					x0, x0 + cell_size[0],
					y0, y0 + cell_size[1],
					z0, z0 + cell_size[2]
				);
			}
		}
	}
}

template<class T, class U>
void Grid<T, U>::set_domain(const geo::Rect<T> & new_domain) {
	domain = new_domain;
}

template<class T, class U>
void Grid<T, U>::insert(const U & node, const Array<T, 3> & pos)
{
	//std::cout << pos[0] << ", " << pos[1] << ", " << pos[2] << std::endl;
	PLB_PRECONDITION(domain.enlarge(min_cell_size).contains_or_on_boundary(pos))
	get_cell(std::floor((pos[0]-x0) / cell_size[0]),
			 std::floor((pos[1]-y0) / cell_size[1]),
			 std::floor((pos[2]-z0) / cell_size[2]))
		.nodes.push_back(node);
}

template<class T, class U>
typename Grid<T, U>::CellType & Grid<T, U>::get_cell(plint i, plint j, plint k)
{
	return cells[k + num_cells.z*(j + num_cells.y*i)];
}

template<class T, class U>
Dot3D Grid<T, U>::get_index(const Array<T, 3> & pos) const
{
	return Dot3D(std::floor((pos[0]-x0) / cell_size[0]),
				 std::floor((pos[1]-y0) / cell_size[1]),
				 std::floor((pos[2]-z0) / cell_size[2]));
}

template<class T, class U>
void Grid<T, U>::clear()
{
	for(typename std::vector<CellType>::iterator it = cells.begin(); it != cells.end(); ++it) {
		it->nodes.clear();
	}
}

} /* namespace fsi */

} /* namespace plb */

#endif /* GRID_HH_ */
