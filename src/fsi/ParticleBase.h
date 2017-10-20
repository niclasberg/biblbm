/*
 * ParticleBase.h
 *
 *  Created on: Jul 22, 2015
 *      Author: niber
 */

#ifndef PARTICLEBASE_H_
#define PARTICLEBASE_H_
#include "core/globalDefs.h"
#include "core/array.h"
#include "core/geometry3D.h"
#include "ParticleShape.h"
#include "ParticleFactory.h"
#include "Quaternion.h"
#include "geometry.h"
#include "Transform.h"
#include "RayTracer.h"

namespace plb {

namespace fsi {

// Forward declarations
template<class T> class Link;
template<class T> class Triangle;
template<class T> class ParticleShapeLibrary;
template<class T> class Boundary;

// Node information object
template<class T>
struct Vertex {
	Array<T, 3> pos;
	Array<T, 3> vel;
	Array<T, 3> force;
};

template<class T>
class ParticleBase3D {
public:
	typedef typename std::vector<Vertex<T> >::iterator vertex_iterator;
	typedef typename std::vector<Vertex<T> >::const_iterator vertex_const_iterator;
	typedef RayTracer<T, Triangle<T>, Vertex<T> > VoxelizerType;

	ParticleBase3D(const ParticleShape<T> *);
	ParticleBase3D(const ParticleBase3D &);
	ParticleBase3D & operator=(const ParticleBase3D &);

	virtual ~ParticleBase3D() { }
	virtual ParticleBase3D * clone() const = 0;

	template<class Derived> bool is_a() const { return this->get_type_id() == Derived::type_id; }
	void copy_nodes_from(const ParticleBase3D &);

	// Getters
	const ParticleShape<T> * shape() const { return shape_; }
	geo::Rect<T> bounding_box() const { return bounding_box_; }
	template<class Arithmetic> geo::Rect<T> bounding_box(const Arithmetic &) const;
	const Array<T, 3> & center_of_mass() const { return center_of_mass_; }
	const T & volume() const { return volume_; }
	const T & area() const { return area_; }

	// Ids
	virtual plint get_type_id() const = 0;
	void set_id(plint id) { this->id = id; }
	plint get_id() const { return id; }
	plint get_shape_id() const { return shape_->get_id(); }

	// Fsi methods
	virtual void reset_forces();
	virtual void compute_forces() = 0;
	virtual void move_vertices() = 0;
	template<class Arithmetic> void map_center_of_mass_to_periodic_grid(const Arithmetic &);

	virtual void update();

	// Setters
	void set_center_of_mass(const Array<T, 3> &);

	// Vertex access
	plint count_nodes() const { return vertices.size(); }
	Vertex<T> & get_node(plint i) { return vertices[i]; }
	const Vertex<T> & get_node(plint i) const { return vertices[i]; }

	// Geometric sampling
	template<class Arithmetic> void get_nodes_in_area(const geo::Rect<T> &, std::vector<Vertex<T> *> &, const Arithmetic &) const;
	void get_nodes_in_area(const geo::Rect<T> &, std::vector<Vertex<T> *> &) const;
	template<class Arithmetic> void append_nodes_in_area(const geo::Rect<T> &, std::vector<Vertex<T> *> &, const Arithmetic &) const;
	void append_nodes_in_area(const geo::Rect<T> &, std::vector<Vertex<T> *> &) const;

	// Voxelizer
	VoxelizerType voxelizer() const { return VoxelizerType(shape_->triangles(), vertices); }

	// Vertex iterators
	vertex_iterator begin() { return vertices.begin(); }
	vertex_const_iterator begin() const { return vertices.begin(); }
	vertex_iterator end() { return vertices.end(); }
	vertex_const_iterator end() const { return vertices.end(); }

	// Serialization for communication and checkpoint
	void pack(std::ostream &) const;
	virtual void pack(std::vector<char> &) const;
	virtual void unpack(char *&);
	virtual void unpack(std::istream &);

	// VTK output
	void write_vtk(FileName) const;
	void write_gmsh(FileName) const;

	// Create shape from configuration
	void store_shape(ParticleShapeLibrary<T> &, std::string tag) const;

	// Lightweight output
	virtual void write_lightweight(std::ostream &) const = 0;

protected:
	T area_, volume_;
	Array<T, 3> center_of_mass_;

private:
	void update_bounding_box();

	const ParticleShape<T> * shape_;
	plint id;
	std::vector<Vertex<T> > vertices;
	geo::Rect<T> bounding_box_;
};

}

}



#endif /* PARTICLEBASE_H_ */
