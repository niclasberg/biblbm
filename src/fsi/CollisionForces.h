/*
 * CollisionForces.h
 *
 *  Created on: Feb 2, 2016
 *      Author: niber
 */

#ifndef COLLISIONFORCES_H_
#define COLLISIONFORCES_H_

namespace plb {

namespace fsi {

namespace collision {

template<class T>
struct CollisionNode {
	Array<T, 3> pos;
	Array<T, 3> * force;
	plint obj_id;
};

/*
 * InteractionFunctional implements the following methods:
 * 	T get_cutoff_distance(): returns the distance after which the force = 0
 * 	T operator()(T distance): returns the force for a given distance
 */
template<class T, class ArithmeticType, class InteractionFunctional>
class CollisionHandler {
public:
	typedef Grid<T, CollisionNode<T> > GridType;

	CollisionHandler(
			const geo::Rect<T> & domain,
			const ArithmeticType & arithmetic,
			const InteractionFunctional & interaction)
	: grid_(domain, 1.),
	  interaction_(interaction),
	  arithmetic_(arithmetic),
	  interaction_distance_sqr_(util::sqr(interaction.get_cutoff_distance()))
	{

	}

	void add_node(const Array<T, 3> & pos, Array<T, 3> * force, plint obj_id)
	{
		CollisionNode<T> node;
		node.force = force;
		node.pos = pos;
		node.obj_id = obj_id;
		grid_.insert(node, pos);
	}

	void compute_collision_forces()
	{
		// Particle-particle interactions
		const Box3D dom = grid_.get_bulk_indices();

		// Do collision for all inner cells
		for(pluint i = dom.x0; i <= dom.x1; ++i) {
			for(pluint j = dom.y0; j <= dom.y1; ++j) {
				for(pluint k = dom.z0; k <= dom.z1; ++k) {
					typename GridType::CellType & cell = grid_.get_cell(i, j, k);

					// Neighbor cells
					typename GridType::CellType * neighbor_cells[7] = {
								&grid_.get_cell(i, j+1, k),
								&grid_.get_cell(i, j+1, k+1),
								&grid_.get_cell(i, j, k+1),
								&grid_.get_cell(i+1, j, k),
								&grid_.get_cell(i+1, j+1, k),
								&grid_.get_cell(i+1, j+1, k+1),
								&grid_.get_cell(i+1, j, k+1)
							};

					for(pluint l = 0; l < cell.nodes.size(); ++l) {
						const CollisionNode<T> & node = cell.nodes[l];

						// Same-cell collisions
						for(pluint m = (l+1); m < cell.nodes.size(); ++m)
							handle_collision(node, cell.nodes[m]);

						// Neighboring cells
						for(pluint n_it = 0; n_it < 7; ++n_it) {
							for(pluint m = 0; m < neighbor_cells[n_it]->nodes.size(); ++m)
								handle_collision(node, neighbor_cells[n_it]->nodes[m]);
						}
					}
				}
			}
		}

		// Treat the case when the domain has periodic edges
		if(arithmetic_.periodic_x()) {
			for(plint j = dom.y0; j <= dom.y1; ++j)
				for(plint k = dom.z0; k <= dom.z1; ++k)
					collide_all_in_cells(grid_.get_cell(dom.x0, j, k), grid_.get_cell(dom.x1, j, k));
		}

		if(arithmetic_.periodic_y()) {
			for(plint i = dom.x0; i <= dom.x1; ++i)
				for(plint k = dom.z0; k <= dom.z1; ++k)
					collide_all_in_cells(grid_.get_cell(i, dom.y0, k), grid_.get_cell(i, dom.y1, k));
		}

		if(arithmetic_.periodic_z()) {
			//std::cout << "Periodic z" << std::endl;
			for(plint i = dom.x0; i <= dom.x1; ++i)
				for(plint j = dom.y0; j <= dom.y1; ++j)
					collide_all_in_cells(grid_.get_cell(i, j, dom.z0), grid_.get_cell(i, j, dom.z1));
		}

		if(arithmetic_.periodic_x() && arithmetic_.periodic_y()) {
			//std::cout << "Periodic xy" << std::endl;
			for(plint k = dom.z0; k <= dom.z1; ++k)
				collide_all_in_cells(grid_.get_cell(dom.x0, dom.y0, k), grid_.get_cell(dom.x1, dom.y1, k));
		}

		if(arithmetic_.periodic_x() && arithmetic_.periodic_z()) {
			//std::cout << "Periodic xz" << std::endl;
			for(plint j = dom.y0; j <= dom.y1; ++j)
				collide_all_in_cells(grid_.get_cell(dom.x0, j, dom.z0), grid_.get_cell(dom.x1, j, dom.z1));
		}

		if(arithmetic_.periodic_y() && arithmetic_.periodic_z()) {
			//std::cout << "Periodic yz" << std::endl;
			for(plint i = dom.x0; i <= dom.x1; ++i)
				collide_all_in_cells(grid_.get_cell(i, dom.y0, dom.z0), grid_.get_cell(i, dom.y1, dom.z1));
		}

		if(arithmetic_.periodic_x() && arithmetic_.periodic_y() && arithmetic_.periodic_z()) {
			//std::cout << "Periodic xyz" << std::endl;
			collide_all_in_cells(grid_.get_cell(dom.x0, dom.y0, dom.z0), grid_.get_cell(dom.x1, dom.y1, dom.z1));
		}
	}

	void compute_wall_collision_forces(const Boundary<T> & boundary)
	{
		const Box3D dom = grid_.get_bulk_indices();

		for(plint iX = dom.x0; iX <= dom.x1; ++iX) {
			for(plint iY = dom.y0; iY <= dom.y1; ++iY) {
				for(plint iZ = dom.z0; iZ <= dom.z1; ++iZ) {
					typename GridType::CellType & cell = grid_.get_cell(iX, iY, iZ);
					if(boundary.does_intersect(cell.domain, interaction_.get_cutoff_distance())) {
						for(plint i = 0; i < cell.nodes.size(); ++i) {
							const CollisionNode<T> & node = cell.nodes[i];

							T dist = boundary.distance_to_boundary(node.pos);
							if(dist <= interaction_.get_cutoff_distance()) {
								*(node.force) -= interaction_(dist) * boundary.get_normal(node.pos);
							}
						}
					}
				}
			}
		}
	}

private:
	void collide_all_in_cells(typename GridType::CellType & cell1, typename GridType::CellType & cell2)
	{
		for(pluint l = 0; l < cell1.nodes.size(); ++l)
			for(pluint m = 0; m < cell2.nodes.size(); ++m)
				handle_collision(cell1.nodes[l], cell2.nodes[m]);
	}

	void handle_collision(const CollisionNode<T> & node, const CollisionNode<T> & node2)
	{
		if(node.obj_id != node2.obj_id) {
			const Array<T, 3> dr = arithmetic_.vec_diff(node2.pos, node.pos);
			const T dist_sqr = normSqr(dr);

			if(dist_sqr <= interaction_distance_sqr_) {
				const T dist =  std::sqrt(dist_sqr);
				const Array<T, 3> df = (-interaction_(dist) / dist) * dr;
				*(node.force) += df;
				*(node2.force) -= df;
			}
		}
	}

private:
	Grid<T, CollisionNode<T> > grid_;
	const InteractionFunctional & interaction_;
	ArithmeticType arithmetic_;
	T interaction_distance_sqr_;
};

} /* namespace collision */

} /* namespace fsi */

} /* namespace plb */



#endif /* COLLISIONFORCES_H_ */
