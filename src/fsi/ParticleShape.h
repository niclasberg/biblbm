#ifndef PARTICLESHAPE_H_
#define PARTICLESHAPE_H_
#include "core/array.h"
#include "Matrix3.h"
#include <vector>
#include "Grid.h"

namespace plb {

namespace fsi {

template<class T> class ParticleShapeLibrary;

/***** ParticleShapeNode *****/
template<class T>
struct Triangle {
	plint i0, i1, i2;
	T area;
	T l_mean;
};

template<class T>
struct Link {
	plint i0, i1;
	plint t1, t2;
	T length;
};

template<class T>
struct AdjacentTrianglePair {
	plint i0, i1, i2, i3;
	T cos_theta0, sin_theta0;
};

struct TriangleQuadruple {
	plint i0, i1, i2, i3, i4, i5;
};

/***** ParticleShape *****/
template<class T>
class ParticleShape {
public:
	typedef typename std::vector<Array<T, 3> >::iterator vertex_iterator;
	typedef typename std::vector<Array<T, 3> >::const_iterator vertex_const_iterator;
	typedef typename std::vector<Triangle<T> >::iterator triangle_iterator;
	typedef typename std::vector<Triangle<T> >::const_iterator triangle_const_iterator;
	typedef typename std::vector<Link<T> >::iterator link_iterator;
	typedef typename std::vector<Link<T> >::const_iterator link_const_iterator;
	typedef typename std::vector<AdjacentTrianglePair<T> >::iterator triangle_pair_iterator;
	typedef typename std::vector<AdjacentTrianglePair<T> >::const_iterator triangle_pair_const_iterator;
	typedef std::vector<TriangleQuadruple>::iterator triangle_quad_iterator;
	typedef std::vector<TriangleQuadruple>::const_iterator triangle_quad_const_iterator;

	~ParticleShape();

	T get_area() const { return area; }
	T get_volume() const { return volume; }
	Matrix<T, 3> get_moment_of_inertia() const { return moment_of_inertia; }
	Array<T, 3> get_center() const { return center; }
	T get_radius() const { return radius; }

	// Identifiers
	pluint get_id() const { return id_; }
	std::string get_tag() const { return tag_; }

	void print_statistics() const;

	// Element counts
	plint count_vertices() const { return vertices_.size(); }
	plint count_triangles() const { return triangles_.size(); }
	plint count_links() const { return links.size(); }

	// Access
	std::vector<Triangle<T> > & triangles() { return triangles_; }
	const std::vector<Triangle<T> > & triangles() const { return triangles_; }

	Triangle<T> & get_triangle(plint i) { return triangles_[i]; }
	const Triangle<T> & get_triangle(plint i) const { return triangles_[i]; }
	Link<T> & get_link(plint i) { return links[i]; }
	const Link<T> & get_link(plint i) const { return links[i]; }
	Array<T, 3> & get_vertex(plint i) { return vertices_[i]; }
	const Array<T, 3> & get_vertex(plint i) const { return vertices_[i]; }

	// Iterators
	vertex_iterator vertices_begin() { return vertices_.begin(); }
	vertex_const_iterator vertices_begin() const { return vertices_.begin(); }
	vertex_iterator vertices_end() { return vertices_.end(); }
	vertex_const_iterator vertices_end() const { return vertices_.end(); }

	triangle_iterator triangles_begin() { return triangles_.begin(); }
	triangle_const_iterator triangles_begin() const { return triangles_.begin(); }
	triangle_iterator triangles_end() { return triangles_.end(); }
	triangle_const_iterator triangles_end() const { return triangles_.end(); }

	link_iterator links_begin() { return links.begin(); }
	link_const_iterator links_begin() const { return links.begin(); }
	link_iterator links_end() { return links.end(); }
	link_const_iterator links_end() const { return links.end(); }

	triangle_pair_iterator triangle_pairs_begin() { return triangle_pairs_.begin(); }
	triangle_pair_iterator triangle_pairs_end() { return triangle_pairs_.end(); }
	triangle_pair_const_iterator triangle_pairs_begin() const { return triangle_pairs_.begin(); }
	triangle_pair_const_iterator triangle_pairs_end() const { return triangle_pairs_.end(); }

	triangle_quad_iterator triangle_quads_begin() { return triangle_quads_.begin(); }
	triangle_quad_iterator triangle_quads_end() { return triangle_quads_.end(); }
	triangle_quad_const_iterator triangle_quads_begin() const { return triangle_quads_.begin(); }
	triangle_quad_const_iterator triangle_quads_end() const { return triangle_quads_.end(); }

private:
	// Make constructor private so only the ShapeLibrary can construct
	ParticleShape(Array<T, 3> *, unsigned int *, unsigned int, unsigned int);

	void compute_connectivity();
	void compute_properties();
	void compute_bounding_volumes();
	void set_id(pluint id) { id_ = id; }
	void set_tag(const std::string & tag) { tag_ = tag; }
	void hash_triangles();
	void find_mesh_neighbors();
	plint hash_position(const Array<T, 3> & pos) const;

	// Underlying mesh
	std::vector<Array<T, 3> > vertices_;
	std::vector<Triangle<T> > triangles_;
	std::vector<Link<T> > links;
	std::vector<AdjacentTrianglePair<T> > triangle_pairs_;
	std::vector<TriangleQuadruple> triangle_quads_;

	// Computed properties
	Matrix<T, 3> moment_of_inertia;		// Area moment of inertia
	T volume;							// Total volume
	T area;
	Array<T, 3> center;					// Volume weighted center
	geo::Rect<T> bounding_box;
	Box3D bounding_box_i;

	// Other stuff
	T radius;							// Bounding sphere radius
	T scale;							// Rescaling coefficient

	// Identifiers
	pluint id_;							// Shape id
	std::string tag_;					// Shape descriptor

	friend class ParticleShapeLibrary<T>;
};

} /* namespace fsi */

} /* namespace plb */

#endif /* PARTICLESHAPE_H_ */
