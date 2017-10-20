/*
 * ImmersedBoundaryDynamics_io.hh
 *
 *  Created on: Apr 20, 2016
 *      Author: niber
 */

#ifndef IMMERSEDBOUNDARYDYNAMICS_IO_HH_
#define IMMERSEDBOUNDARYDYNAMICS_IO_HH_
#include "ImmersedBoundaryDynamics.h"

/****** Io ******/
namespace plb {

namespace fsi {
template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::save_checkpoint(FileName fname) const
{
	// Header file name
	FileName header_file = fname;
	std::string base_folder = global::directories().getImageOutDir();
	if( ! base_folder.empty() && base_folder.at(base_folder.length()-1) != '/')
		base_folder.append("/");
	header_file.setPath(base_folder);
	header_file.setExt("plb");

	// Processor file base path
	FileName proc_file;
	proc_file.setPath(header_file.getName() + "/");
	proc_file.setExt("plb");

	// Generate list of file names for all procs
	std::vector<FileName> file_names;
	for(plint i = 0; i < global::mpi().getSize(); ++i) {
		FileName fn = proc_file;
		fn.setName(util::val2str(i));
		file_names.push_back(fn);
	}

	// Write header
	if(global::mpi().isMainProcessor()) {
		std::ofstream header_out(header_file.get().c_str());

		if(! header_out) {
			std::cerr << "Could not open the checkpoint header file " << header_file.get() << " for writing" << std::endl;
			exit(-1);
		}

		plint num_procs = global::mpi().getSize();
		header_out << (num_particles) << std::endl;
		header_out << num_procs << std::endl;

		for(plint i = 0; i < file_names.size(); ++i)
			header_out << file_names[i].get() << std::endl;
	}

	FileName fname_local = file_names[global::mpi().getRank()];
	fname_local.setPath(base_folder + fname_local.getPath());

	io::mkdir(fname_local.getPath().c_str());

	std::ofstream out(fname_local.get().c_str(), std::ios::binary);
	if(! out) {
		std::cerr << "Could not open the file " << fname_local.get() << " for writing" << std::endl;
		exit(-1);
	}

	utils::pack<plint>(out, particles.size());
	for(typename ObjMapType::const_iterator it = particles.begin(); it != particles.end(); ++it)
		it->second->pack(out);
	out.close();
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::load_checkpoint(FileName fname)
{
	// Clear particles and nodes
	clear();

	// Header file
	FileName header_file = fname;
	std::string base_folder = global::directories().getInputDir();
	if( ! base_folder.empty() && base_folder.at(base_folder.length()-1) != '/')
		base_folder.append("/");
	header_file.setPath(base_folder);
	header_file.setExt("plb");

	// Read header file
	plint num_particles, num_files;
	std::ifstream * in;
	std::vector<std::string> filenames;
	if(global::mpi().isMainProcessor()) {
		std::ifstream header(header_file.get().c_str());
		if( ! header) {
			std::cerr << "Could not open the checkpoint header file " << header_file.get() << " for reading" << std::endl;
			exit(-1);
		}

		header >> num_particles;
		header >> num_files;

		std::string file_name;
		std::getline(header, file_name);
		for(plint i = 0; i < num_files; ++i) {
			std::getline(header, file_name);

			filenames.push_back(base_folder + file_name);
		}

		header.close();
	}

	// Broadcast the total number of particles and the number of files read
	global::mpi().bCast(&num_particles, 1);
	global::mpi().bCast(&num_files, 1);

	this->num_particles = num_particles;

	// Read files and broadcast the content
	std::vector<char> buffer;
	ParticleType * particle;
	plint buffer_size;
	plint particles_read;

	for(plint i = 0; i < num_files; ++i) {
		buffer.clear();
		if(global::mpi().isMainProcessor()) {
			std::ifstream file(filenames[i].c_str(), std::ios::binary);
			if( ! file.good()) {
				std::cerr << "Could not open the file " << filenames[i] << " for reading" << std::endl;
				exit(-1);
			}

			buffer.insert(buffer.begin(), std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>());
			buffer_size = buffer.size();
			file.close();
		}

		global::mpi().bCast(&buffer_size, 1);

		if( ! global::mpi().isMainProcessor()) {
			buffer.resize(buffer_size);
		}

		global::mpi().bCast(&(buffer[0]), buffer_size);

		char * it = &(buffer[0]);
		utils::unpack(it, particles_read);

		// Create particles
		for(plint j = 0; j < particles_read; ++j) {
			particle = particleFactory<T>().create(it, shape_library);
			add_particle(particle, particle->get_id());
			delete particle;
		}
	}
	is_init_ = false;
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::write_particles_as_vtk(pluint iter)
{
	for(ObjMapIterator it = particles_begin(); it != particles_end(); ++it) {
		// Create file name
		std::stringstream ss;
		ss << "obj" << it->first << "-" << iter;

		FileName fname;
		fname.setName(ss.str());
		fname.setPath(global::directories().getOutputDir());
		fname.setExt(".vtu");

		it->second->write_vtk(fname.get());
	}
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::write_lightweight_particle_data(pluint it)
{
	// Create output file
	std::stringstream ss;
	ss << "objsP" << global::mpi().getRank() << "-" << it;;

	FileName fname;
	fname.setName(ss.str());
	fname.setPath(global::directories().getOutputDir());
	fname.setExt(".txt");

	std::ofstream out(fname.get().c_str(), std::ios::out);
	out << "iteration:" << it << std::endl;

	for(ObjMapConstIterator it = particles_begin(); it != particles_end(); ++it) {
		ParticleBase3D<T> * p = it->second;
		p->write_lightweight(out);
		out << std::endl;
	}

	out.close();
}

} /* namespace fsi */
} /* namespace plb */

#endif /* IMMERSEDBOUNDARYDYNAMICS_IO_HH_ */
