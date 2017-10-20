/*
 * ParticleBase.hh
 *
 *  Created on: Jul 22, 2015
 *      Author: niber
 */

#ifndef PARTICLEBASE_HH_
#define PARTICLEBASE_HH_
#include "ParticleBase.h"
#include "utils.h"
#include "io/base64.h"
#include "ParticleShapeFactory.h"

namespace plb {

namespace fsi {

/*** Constructors ***/
template<class T>
ParticleBase3D<T>::ParticleBase3D(
		const ParticleShape<T> * shape)
: shape_(shape),
  id(-1),
  area_(shape->get_area()),
  volume_(shape->get_volume()),
  center_of_mass_(shape->get_center())
{
	// Create vertices
	vertices.resize(this->shape()->count_vertices());

	plint i = 0;
	for(typename ParticleShape<T>::vertex_const_iterator it = this->shape()->vertices_begin();
			it != this->shape()->vertices_end(); ++it) {
		vertices[i].pos = *it;
		vertices[i].vel.resetToZero();
		vertices[i].force.resetToZero();
		++i;
	}

	update_bounding_box();
}

template<class T>
void ParticleBase3D<T>::update()
{
	update_bounding_box();
}

template<class T>
void ParticleBase3D<T>::update_bounding_box()
{
	bounding_box_.x0 = std::numeric_limits<T>::max();
	bounding_box_.y0 = std::numeric_limits<T>::max();
	bounding_box_.z0 = std::numeric_limits<T>::max();
	bounding_box_.x1 = -std::numeric_limits<T>::max();
	bounding_box_.y1 = -std::numeric_limits<T>::max();
	bounding_box_.z1 = -std::numeric_limits<T>::max();

	for(vertex_const_iterator it = begin(); it != end(); ++it) {
		bounding_box_.x0 = std::min(bounding_box_.x0, it->pos[0]);
		bounding_box_.y0 = std::min(bounding_box_.y0, it->pos[1]);
		bounding_box_.z0 = std::min(bounding_box_.z0, it->pos[2]);
		bounding_box_.x1 = std::max(bounding_box_.x1, it->pos[0]);
		bounding_box_.y1 = std::max(bounding_box_.y1, it->pos[1]);
		bounding_box_.z1 = std::max(bounding_box_.z1, it->pos[2]);
	}
}

template<class T>
ParticleBase3D<T>::ParticleBase3D(const ParticleBase3D & rhs)
: shape_(rhs.shape_),
  vertices(rhs.vertices),
  id(rhs.id),
  //voxelizer_(shape_->triangles(), vertices),
  area_(rhs.area_),
  volume_(rhs.volume_),
  center_of_mass_(rhs.center_of_mass_)
{
	update();
}

template<class T>
ParticleBase3D<T> & ParticleBase3D<T>::operator=(const ParticleBase3D & rhs)
{
	if(&rhs != this) {
		shape_ = rhs.shape_;
		vertices = rhs.vertices;
		id = rhs.id;
		area_ = rhs.area_;
		volume_ = rhs.volume_;
		center_of_mass_ = rhs.center_of_mass_;
		update();
	}
	return *this;
}

template<class T>
void ParticleBase3D<T>::copy_nodes_from(const ParticleBase3D & rhs)
{
	assert(vertices.size() == this->shape()->count_vertices());
	vertices = rhs.vertices;
	update();
}


/*** Vertex transformations and access ***/
template<class T>
	template<class Arithmetic>
void ParticleBase3D<T>::map_center_of_mass_to_periodic_grid(const Arithmetic & arithmetic)
{
	Array<T, 3> pos = center_of_mass_;
	arithmetic.remap_position(pos);
	set_center_of_mass(pos);
}

template<class T>
void ParticleBase3D<T>::set_center_of_mass(const Array<T, 3> & pos)
{
	const Array<T, 3> delta = pos - center_of_mass_;

	// Adjust the position of the vertices
	for(vertex_iterator it = begin(); it != end(); ++it)
		it->pos += delta;

	center_of_mass_ = pos;

	// Adjust bounding box
	bounding_box_.x0 += delta[0];
	bounding_box_.x1 += delta[0];
	bounding_box_.y0 += delta[1];
	bounding_box_.y1 += delta[1];
	bounding_box_.z0 += delta[2];
	bounding_box_.z1 += delta[2];
}

template<class T>
void ParticleBase3D<T>::get_nodes_in_area(const geo::Rect<T> & area, std::vector<Vertex<T> *> & ret) const
{
	ret.clear();
	append_nodes_in_area(area, ret);
}

template<class T>
	template<class Arithmetic>
void ParticleBase3D<T>::append_nodes_in_area(const geo::Rect<T> & area, std::vector<Vertex<T> *> & ret, const Arithmetic & arithmetic) const
{
	for(vertex_const_iterator it = begin(); it != end(); ++it) {
		if(area.contains(it->pos, arithmetic))
			ret.push_back(const_cast<Vertex<T> *>(&(*it)));
	}
}

template<class T>
	template<class Arithmetic>
void ParticleBase3D<T>::get_nodes_in_area(const geo::Rect<T> & area, std::vector<Vertex<T> *> & ret, const Arithmetic & arithmetic) const
{
	ret.clear();
	append_nodes_in_area(area, ret, arithmetic);
}

template<class T>
void ParticleBase3D<T>::append_nodes_in_area(const geo::Rect<T> & area, std::vector<Vertex<T> *> & ret) const
{
	for(vertex_const_iterator it = begin(); it != end(); ++it) {
		if(area.contains(it->pos))
			ret.push_back(const_cast<Vertex<T> *>(&(*it)));
	}
}

/*********** Utility *******/
template<class T>
	template<class Arithmetic>
geo::Rect<T> ParticleBase3D<T>::bounding_box(const Arithmetic & arithmetic) const
{
	geo::Rect<T> ret = bounding_box_;
	arithmetic.remap_position_x(ret.x0);
	arithmetic.remap_position_x(ret.x1);
	arithmetic.remap_position_y(ret.y0);
	arithmetic.remap_position_y(ret.y1);
	arithmetic.remap_position_z(ret.z0);
	arithmetic.remap_position_z(ret.z1);
}

/*********** Fsi ***********/
template<class T>
void ParticleBase3D<T>::reset_forces()
{
	// Reset forces to zero
	for(vertex_iterator it = this->begin(); it != this->end(); ++it)
		it->force.resetToZero();
}

/****** Convert to shape **/
template<class T>
void ParticleBase3D<T>::store_shape(ParticleShapeLibrary<T> & shape_library, std::string tag) const
{
	std::vector<Array<T, 3> > vertices;
	std::vector<unsigned int> indices;
	for(vertex_const_iterator it = begin(); it != end(); ++it)
		vertices.push_back(it->pos);
	for(typename ParticleShape<T>::triangle_const_iterator it = shape()->triangles_begin();
			it != shape()->triangles_end(); ++it) {
		indices.push_back(it->i0);
		indices.push_back(it->i1);
		indices.push_back(it->i2);
	}

	shape_library.store_mesh(&(vertices[0]), &(indices[0]), vertices.size(), shape_->count_triangles(), tag);
}


/*********** Io ***********/
template<class T>
void ParticleBase3D<T>::pack(std::vector<char> & buff) const
{
	utils::pack(buff, this->get_type_id());
	utils::pack(buff, this->shape()->get_id());
	utils::pack(buff, this->get_id());
}

template<class T>
void ParticleBase3D<T>::pack(std::ostream & out) const
{
	std::vector<char> v;
	pack(v);
	out.write(&(v[0]), v.size());
}

template<class T>
void ParticleBase3D<T>::unpack(char *& buff)
{
	plint id;
	utils::unpack(buff, id);
	this->set_id(id);
}

template<class T>
void ParticleBase3D<T>::unpack(std::istream & in)
{
	plint id;
	utils::unpack(in, id);
	this->set_id(id);
}


namespace detail {
	template<class T>
	class Writer {
	public:
		Writer(std::ostream & out_) : out(out_) { }
		template<class T2>
		void convert_and_write(T2 val)
		{
			T v2 = val;
			out.write(reinterpret_cast<char *>(&v2), sizeof(T));
		}
	private:
		std::ostream & out;
	};
}

template<class T>
void ParticleBase3D<T>::write_vtk(FileName fname) const
{
	std::string uint_name, int_name, float_name;
	switch(sizeof(float)) {
	case 4:
		float_name = "Float32";
		break;
	case 8:
		float_name = "Float64";
		break;
	default:
		float_name = "Unsupported";
		break;
	}

	switch(sizeof(unsigned int)) {
	case 4:
		uint_name = "UInt32";
		break;
	case 8:
		uint_name = "UInt64";
		break;
	default:
		uint_name = "Unsupported";
		break;
	}

	switch(sizeof(int)) {
	case 4:
		int_name = "Int32";
		break;
	case 8:
		int_name = "Int64";
		break;
	default:
		uint_name = "Unsupported";
		break;
	}

	std::ofstream out(fname.get().c_str(), std::ios::out|std::ios::trunc|std::ios::binary);

	if(!out.good()) {
		std::cerr << "Warning: DeformableParticle<T>::write_vtk: Could not open the file "
			<< fname.get() << " for writing" << std::endl;
	}

	// Pack position data
	pluint num_stuff = 5;
	pluint sizes[5] = {
				3 * vertices.size() * sizeof(float), 				// positions
				3 * vertices.size() * sizeof(float), 				// forces
				3 * this->shape()->count_triangles() * sizeof(int),	// connectivity
				this->shape()->count_triangles() * sizeof(int), 	// offsets
				this->shape()->count_triangles() * sizeof(int)		// types
	};
	pluint offsets[5];

	// compute offsets
	offsets[0] = 0;
	for(plint i = 1; i < num_stuff; ++i)
		offsets[i] = offsets[i-1] + sizeof(unsigned int) + sizes[i-1];

	out << "<?xml version=\"1.0\"?>\n";
	out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"" << uint_name << "\" >\n";
	out << "  <UnstructuredGrid>\n";
	out << "    <Piece NumberOfPoints=\""<< this->vertices.size() <<"\" NumberOfCells=\""<< shape()->count_triangles() << "\" >" << std::endl;
	out << "      <Points>" << std::endl;
	out << "        <DataArray type=\""<< float_name << "\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offsets[0] << "\" />" << std::endl;
	out << "      </Points>" << std::endl;
	out << "      <PointData>" << std::endl;
	out << "        <DataArray type=\"" << float_name << "\" NumberOfComponents=\"3\" Name=\"force_tot\" format=\"appended\" offset=\"" << offsets[1] << "\" />" << std::endl;
	out << "      </PointData>" << std::endl;
	out << "      <Cells>" << std::endl;
	out << "        <DataArray type=\"" << int_name << "\" Name=\"connectivity\" format=\"appended\" offset=\"" << offsets[2] << "\" />" << std::endl;
	out << "        <DataArray type=\"" << int_name << "\" Name=\"offsets\" format=\"appended\" offset=\"" << offsets[3] << "\" />" << std::endl;
	out << "        <DataArray type=\"" << int_name << "\" Name=\"types\" format=\"appended\" offset=\"" << offsets[4] << "\" />" << std::endl;
	out << "      </Cells>" << std::endl;
	out << "    </Piece>" << std::endl;
	out << "  </UnstructuredGrid>" << std::endl;

	// Write appended data
	out << "  <AppendedData encoding=\"raw\">" << std::endl;
	out << "_";

	// Writers (convenience wrappers around the stream object taking care of variable conversions)
	detail::Writer<float> float_writer(out);
	detail::Writer<unsigned int> uint_writer(out);
	detail::Writer<int> int_writer(out);

	// Vertices
	uint_writer.convert_and_write(sizes[0]);
	for(plint i = 0; i < vertices.size(); ++i) {
		float_writer.convert_and_write(vertices[i].pos[0]);
		float_writer.convert_and_write(vertices[i].pos[1]);
		float_writer.convert_and_write(vertices[i].pos[2]);
	}

	// force data
	uint_writer.convert_and_write(sizes[1]);
	for(plint i = 0; i < vertices.size(); ++i) {
		float_writer.convert_and_write(vertices[i].force[0]);
		float_writer.convert_and_write(vertices[i].force[1]);
		float_writer.convert_and_write(vertices[i].force[2]);
	}

	// Connectivity data
	uint_writer.convert_and_write(sizes[2]);
	for(plint i = 0; i < shape()->count_triangles(); ++i) {
		int_writer.convert_and_write(shape()->get_triangle(i).i0);
		int_writer.convert_and_write(shape()->get_triangle(i).i1);
		int_writer.convert_and_write(shape()->get_triangle(i).i2);
	}

	// Offsets
	uint_writer.convert_and_write(sizes[3]);
	for(plint i=0; i < shape()->count_triangles(); i++)
		int_writer.convert_and_write((i+1)*3);

	// Pack types
	uint_writer.convert_and_write(sizes[4]);
	for(plint i=0; i < shape()->count_triangles(); i++)
		int_writer.convert_and_write(5);

	out << std::endl << "  </AppendedData>" << std::endl;
	out << "</VTKFile>" << std::endl;

	out.close();
}

template<class T>
void ParticleBase3D<T>::write_gmsh(FileName fname) const
{
	char version[] = "2.1";
	int file_type = 0;
	int data_size = sizeof(double);
	int element_type = 2; // Triangles
	std::vector<int> tags;
	tags.push_back(0); tags.push_back(1); tags.push_back(0);

	std::ofstream out(fname.get().c_str(), std::ios::out|std::ios::trunc);
	out << "$MeshFormat" << std::endl;
	// Version file_type data_size
	out << version << " " << file_type << " " << data_size << std::endl;
	out << "$EndMeshFormat" << std::endl;

	// Write nodes
	out << "$Nodes" << std::endl;
	out << count_nodes() << std::endl;
	for(int i = 0; i < count_nodes(); ++i)
		out << i << " " << get_node(i).pos[0] << " " << get_node(i).pos[1] << " " << get_node(i).pos[2] << std::endl;
	out << "$EndNodes" << std::endl;

	// Write elements
	out << "$Elements" << std::endl;
	out << shape()->count_triangles() << std::endl;

	for(int i = 0; i < shape()->count_triangles(); ++i) {
		out << i << " " << 2 << " " << 3 << " ";
		for(int j = 0; j < tags.size(); ++j)
			out << tags[j] << " ";
		out << shape_->get_triangle(i).i0 << " " << shape_->get_triangle(i).i1 << " " << shape_->get_triangle(i).i2;
		out << std::endl;
	}
	out << "$EndElements";

	out.close();
}

}

}



#endif /* PARTICLEBASE_HH_ */
