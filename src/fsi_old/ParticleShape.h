#ifndef PARTICLESHAPE_H_
#define PARTICLESHAPE_H_
#include "core/array.h"
#include "Matrix3.h"
#include <vector>
#include "Grid.h"

namespace plb {

namespace fsi {

template<class T> class ParticleShapeLibrary;

/***** Particle shape descriptors *****/
/*template<class T>
class ParticleShapeDescriptor {
public:
	virtual bool point_is_inside(const Array<T, 3> &) const
	{
		PLB_ASSERT(false);
	}
};

template<class T>
class EllipsoidalShapeDescriptor : public ParticleShapeDescriptor<T> {
public:
	virtual bool point_is_inside(const Array<T, 3> & p) const
	{
		return util::sqr(p[0]/a) + util::sqr(p[1]/b) + util::sqr(p[2]/c) < 1;
	}
private:
	T a, b, c;
};

template<class T>
class BiconcaveShapeDescriptor : public ParticleShapeDescriptor<T> {
public:
	virtual bool point_is_inside(const Array<T, 3> & p) const
	{
		const T rho_sqr = util::sqr(p[0]) + util::sqr(p[1]);
		const T z_sqr = util::sqr(p[2]);
		return util::sqr(rho_sqr + z_sqr) + P*rho_sqr + Q*z_sqr + R < 0;
	}
private:
	T a, b, c;
};*/

/***** ParticleShapeNode *****/
template<class T>
struct ParticleShapeNode {
	Array<T, 3> centroid;			// Centroid
	Array<T, 3> vertices[3];		// Vertices in the triangle
	Array<T, 3> normal;				// Normal vector
	T lambda;						// Curvature estimate
	T area;							// Surface area
	std::vector<plint> neighbors;	// Neighboring node indices
	plint id;
};

/***** ParticleShape *****/
template<class T>
class ParticleShape {
public:
	ParticleShape(Array<T, 3> *, unsigned int *, unsigned int, unsigned int);
	~ParticleShape();

	const Array<T, 3> * get_vertex_ptr() const { return vertices_; }
	const unsigned int * get_index_ptr() const { return indices; }
	unsigned int get_num_elements() const { return num_triangles; }
	unsigned int get_num_vertices() const { return num_vertices; }

	ParticleShapeNode<T> & get_node(pluint i) { return nodes[i]; }
	const ParticleShapeNode<T> & get_node(pluint i) const { return nodes[i]; }

	Array<T, 3> & get_centroid(pluint i) { return nodes[i].centroid; }
	const Array<T, 3> & get_centroid(pluint i) const { return nodes[i].centroid; }

	T get_volume() const { return volume; }
	Matrix<T, 3> get_moment_of_inertia() const { return moment_of_inertia; }
	Array<T, 3> get_center() const { return center; }

	T get_radius() const { return radius; }
	pluint get_id() const { return id; }

	plint get_closest_element_facing_in_direction(const Array<T, 3> &, const Array<T, 3> &, T cutoff = std::numeric_limits<T>::max()) const;

	Array<T, 3> get_average_distance_within_sphere(const Array<T, 3> &, T) const;

	void trace_ray(const Array<T, 3> & pos, const Array<T, 3> & dir, plint &, T &) const;
	T get_closest_distance_to_surface(const Array<T, 3> &) const;
	T get_closest_distance_to_surface(const Array<T, 3> &, const Array<T, 3> &) const;

	//bool does_overlap(const ParticleShape &, const Array<T, 3> &, const Quaternion<T> &) const;

private:
	void compute_properties();
	void compute_bounding_volumes();
	void set_id(pluint id_) { id = id_; }
	void hash_triangles();
	void find_mesh_neighbors();
	plint hash_position(const Array<T, 3> & pos) const;

	// Underlying mesh
	unsigned int num_triangles, num_vertices;
	Array<T, 3> * vertices_;
	unsigned int * indices;

	// Nodes
	std::vector<ParticleShapeNode<T> > nodes;

	// Grid with the set of points in the shape (for accelerated collision detection)
	std::vector<std::vector<ParticleShapeNode<T> *> > triangle_hash;
	Dot3D hash_dim;

	// Computed properties
	Matrix<T, 3> moment_of_inertia;		// Area moment of inertia
	T volume;							// Total volume
	Array<T, 3> center;					// Volume weighted center
	geo::Rect<T> bounding_box;
	Box3D bounding_box_i;

	// Other stuff
	T radius;							// Bounding sphere radius
	pluint id;							// Shape id
	T scale;							// Rescaling coefficient

	friend class ParticleShapeLibrary<T>;
};

} /* namespace fsi */

} /* namespace plb */

#endif /* PARTICLESHAPE_H_ */
