/*
 * ParticleShapeFactory.hh
 *
 *  Created on: Mar 4, 2014
 *      Author: niber
 */

#ifndef PARTICLESHAPEFACTORY_HH_
#define PARTICLESHAPEFACTORY_HH_

namespace plb {

namespace fsi {

namespace detail {

int ReadToPos(FILE *fp, std::string str)
{
  char line[512];
  while (!feof(fp)) {
      fgets(line, 512, fp);
      if (std::strstr(line, str.c_str()) != NULL) {
          return 1;
      }
  }
  return 0;
}

template<class T>
void read_gmsh_mesh(std::string fname, unsigned int **nodpek, Array<T, 3> ** coord, unsigned int *elnum, unsigned int *nnode, unsigned int *nobf);

void read_gmsh_mesh(std::string fname, unsigned int **nodpek, Array<double,3> ** coord, unsigned int *elnum, unsigned int *nnode, unsigned int *nobf)
{
	unsigned int tnnode, telnum, ngrp, nseg;
	char ctmp[256];
	Array<double, 3> * fcoord;
	unsigned int *fnodpek;
	unsigned int i, j, k, pr, n1, n2, n3, nn, mode, un, er;
	double x, y, z;

	FILE *fp = fopen(fname.c_str(), "r");
	if (fp == NULL) {
		std::cerr << "ERROR: the input file \"" << fname << "\" not found\n";
		exit(-1);
	} else {
		if (ReadToPos(fp, "$Nodes")) {
			fscanf(fp, "%i", &tnnode);

			fcoord = new Array<double, 3>[tnnode];

			for (i = 0; i < tnnode; i++) {
				fscanf(fp, "%d%lf%lf%lf", &j, &x, &y, &z);
				j--;
				fcoord[j][0] = x;
				fcoord[j][1] = y;
				fcoord[j][2] = z;
			}

			fgets(ctmp, 256, fp);
			fgets(ctmp, 256, fp);
			fgets(ctmp, 256, fp); // $Elements
			fscanf(fp, "%d", &telnum);

			fnodpek = new unsigned int[telnum * 3];

			for (i = 0; i < telnum; i++) {
				fscanf(fp, "%i%i%i%i%i%i%i%i%i", &j, &mode, &nn, &pr, &er, &un,
						&n1, &n2, &n3);
				j--;
				fnodpek[j * 3] = n1 - 1;
				fnodpek[j * 3 + 1] = n2 - 1;
				fnodpek[j * 3 + 2] = n3 - 1;
			}
		}
		fclose(fp);

		*nodpek = fnodpek;
		*coord = fcoord;
		*nobf = 3;
		*elnum = telnum;
		*nnode = tnnode;
	}
}

} /* namespace detail */


template<class T>
ParticleShapeLibrary<T>::~ParticleShapeLibrary()
{
	clear();
}

template<class T>
void ParticleShapeLibrary<T>::clear()
{
	// Free the allocated memory
	for(typename std::vector<ParticleShape<T> *>::iterator it = id_shape_map.begin(); it != id_shape_map.end(); ++it)
		delete *it;
	id_shape_map.clear();
	name_id_map.clear();
}

template<class T>
void ParticleShapeLibrary<T>::read_and_store_mesh(const std::string & file_name, const std::string & tag)
{
	unsigned int num_els, num_nodes, nodes_per_el;
	Array<T, 3> * vertices;
	unsigned int * indices;

	detail::read_gmsh_mesh(file_name, &indices, &vertices, &num_els, &num_nodes, &nodes_per_el);
	store_mesh(vertices, indices, num_nodes, num_els, tag);

	// free memory
	delete [] indices;
	delete [] vertices;
}

template<class T>
void ParticleShapeLibrary<T>::store_mesh(
		Array<T, 3> * vertices,
		unsigned int * indices,
		unsigned int num_vertices,
		unsigned int num_cells,
		const std::string & tag)
{
	PLB_PRECONDITION(name_id_map.find(tag) == name_id_map.end())

	ParticleShape<T> * shape = new ParticleShape<T>(vertices, indices, num_vertices, num_cells);
	plint next_id = id_shape_map.size();
	shape->set_id(next_id);
	shape->set_tag(tag);

	// Store by tag and id
	id_shape_map.push_back(shape);
	name_id_map[tag] = next_id;
}

template<class T>
const ParticleShape<T> * ParticleShapeLibrary<T>::get_by_id(plint id) const
{
	//std::cout << id << " / " << id_shape_map.size() << std::endl;
	assert(id < id_shape_map.size());
	return id_shape_map[id];
}

template<class T>
ParticleShape<T> * ParticleShapeLibrary<T>::get_by_id(plint id)
{
	assert(id_shape_map.find(id) != id_shape_map.end());
	return id_shape_map[id];
}

template<class T>
const ParticleShape<T> * ParticleShapeLibrary<T>::get_by_tag(std::string tag) const
{
	return id_shape_map[get_id_by_tag(tag)];
}

template<class T>
plint ParticleShapeLibrary<T>::get_id_by_tag(std::string tag) const
{
	if(name_id_map.find(tag) == name_id_map.end()) {
		std::cerr << "Could not find the ParticleShape " << tag << std::endl;
		exit(-1);
	}
	return name_id_map.at(tag);
}

template<class T>
ParticleShape<T> * ParticleShapeLibrary<T>::get_by_tag(std::string tag)
{
	return id_shape_map[get_id_by_tag(tag)];
}

template<class T>
T ParticleShapeLibrary<T>::get_max_particle_radius() const
{
	T ret = 0;
	for(typename std::vector<ParticleShape<T> *>::const_iterator it = id_shape_map.begin(); it != id_shape_map.end(); ++it)
		ret = std::max(ret, (*it)->get_radius());
	return ret;
}

template<class T>
plint ParticleShapeLibrary<T>::get_shape_count() const
{
	return id_shape_map.size();
}

} /* namespace fsi */

} /* namespace plb */


#endif /* PARTICLESHAPEFACTORY_HH_ */
