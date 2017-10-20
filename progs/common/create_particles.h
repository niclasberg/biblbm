#ifndef CREATE_PARTICLES_H_
#define CREATE_PARTICLES_H_
#include "fsi/headers.h"

// Create an initial particle distribution by initializing the
// particles randomly at a fraction of their size. The particles
// are then grown to their full size while under the influence of
// repulsive forces.
template<class T>
void create_particles(
		ImmersedBoundaryDynamics3D<T, DESCRIPTOR, Periodicity> & fsi,
		MultiBlockLattice3D<T,DESCRIPTOR> & lattice)
{
	// Create intial condition for particles
	if(load_checkpoint) {
		pcout << "Loading particle checkpoint" << std::endl;
		fsi.load_checkpoint("checkpoint_fsi");
	} else {
		plint Nit = 40000;
		T initial_scale = 0.2;

		pcout << "Generating initial distribution" << std::endl;

		// Create a new immersed boundary dynamics object for the rigid body simulation
		ParticleShapeLibrary<T> shape_library2;
		ImmersedBoundaryDynamics3D<T, DESCRIPTOR, Periodicity> fsi2(lattice, shape_library2);

		// Create a boundary that is a little bit smaller than the actual one,
		// this ensures that all particles are within the domain
		PipeBoundary<T> boundary2 = pipe_boundary;
		boundary2.radius = pipe_boundary.radius - 1.5;

		plint num_particles = 0;
		T max_radius = 0;
		std::vector<RigidParticle3D<T> *> rigidParticles;
		std::vector<T> particleInitialRadii;
		for(int i = 0; i < particleTypes.size(); ++i) {
			// Compute the total number of particles
			num_particles += particleTypes[i].numParticles;

			// Store the equilibrium shape
			particleTypes[i].particle->store_shape(shape_library2, particleTypes[i].name.c_str());
			ParticleShape<T> * rigid_shape = shape_library2.get_by_tag(particleTypes[i].name.c_str());

			// Create a rigid particle with the same shape but 90% smaller
			RigidParticle3D<T> * rigid_particle = new RigidParticle3D<T>(rigid_shape);
			rigid_particle->scale() = initial_scale;
			rigid_particle->density() = 1.;
			rigid_particle->damping() = 0.1;
			rigid_particle->update();
			rigidParticles.push_back(rigid_particle);

			// Update max radius
			T initialRadius = initial_scale*rigid_shape->get_radius() + 1.;
			particleInitialRadii.push_back(initialRadius);
			max_radius = std::max(max_radius, initialRadius);
		}

		// Generate initial distribution
		std::vector<Array<T, 3> > positions;
		std::vector<Quaternion<T> > orientations;
		if(global::mpi().isMainProcessor()) {
			// Generate positions
			ParticlePositionInitializer<T> position_initializer(
					lattice.getBoundingBox(),
					&boundary2, max_radius);

			for(int i = 0; i < particleTypes.size(); ++i)
				position_initializer.generate_points(particleTypes[i].numParticles, particleInitialRadii[i], i);

			// Get positions
			for(int i = 0; i < particleTypes.size(); ++i) {
				std::vector<Array<T, 3> > posTmp;
				position_initializer.get_points(i, posTmp);
				positions.insert(positions.end(), posTmp.begin(), posTmp.end());
			}

			// Generate random orientations
			for(plint i = 0; i < num_particles; ++i) {
				Array<T, 3> dir((T)std::rand() / RAND_MAX,
						(T)std::rand() / RAND_MAX,
						(T)std::rand() / RAND_MAX);
				T dir_norm = norm(dir);
				if(dir_norm > std::numeric_limits<T>::epsilon()) {
					T theta = 2. * M_PI * (T)std::rand() / RAND_MAX;
					orientations.push_back(fsi::Quaternion<T>(theta, dir/dir_norm));
				} else {
					orientations.push_back(fsi::Quaternion<T>());
				}
			}
		} else {
			positions.resize(num_particles);
			orientations.resize(num_particles);
		}

		// Broadcast to all processors
		global::mpi().bCast(&initial_scale, 1);
		global::mpi().bCast(reinterpret_cast<T*>(&(positions[0])), num_particles*3);
		global::mpi().bCast(reinterpret_cast<T*>(&(orientations[0])), num_particles*4);

		pcout << "Initial distribution created with the RBCs at " << 100*initial_scale << "% of their full size" << std::endl;
		pcout << "Growing particles to their full size" << std::endl;

		// Create and insert particles to the fsi object
		{
			int N;
			for(int i = 0, N = 0; i < particleTypes.size(); ++i)
				for(plint j = 0; j < particleTypes[i].numParticles; ++j, ++N) {
					RigidParticle3D<T> * p = rigidParticles[i];
					p->scale() = initial_scale;
					p->set_center_of_mass(positions[N]);
					p->orientation() = orientations[N];
					p->update();
					fsi2.add_particle(p);
				}
		}

		// Initialize fsi
		fsi2.init();

		// Potential force
		SpringPotential<T> potential(1., 5e-2);

		// Grow the particles while under the influence of collision forces
		T scale;
		pipe_boundary.writeVTK(global::directories().getOutputDir() + std::string("boundary.vtu"));
		for(plint it = 0; it <= Nit; ++it) {
			if((it % 100) == 0)
				pcout << "  Iteration " << it << " / " << Nit << std::endl;
			//if((it % 500) == 0) fsi2.write_particles_as_vtk(it);

			if((it % 10) == 0) {
				// Rescale all particles
				T fraction = (T)it / (T) Nit;
				scale = initial_scale + (1. - initial_scale)*fraction;

				for(typename ImmersedBoundaryDynamics3D<T, DESCRIPTOR, Periodicity>::ObjMapIterator iter = fsi2.particles_begin();
						iter != fsi2.particles_end(); ++iter) {
					plb::fsi::RigidParticle3D<T> * p = dynamic_cast<plb::fsi::RigidParticle3D<T> *>(iter->second);
					if(p) {
						p->scale() = scale;
						p->update();
					}
				}

				fsi2.synchronize_particle_states();
			}

			fsi2.compute_collision_forces(potential, &boundary2);
			fsi2.move_vertices();
		}

		// Scatter all particle positions and orientations
		CommunicationBuffer comm_buff;
		std::vector<plint> proc_list;
		for(plint i = 0; i < global::mpi().getSize(); ++i)
			if(i != global::mpi().getRank())
				proc_list.push_back(i);
		comm_buff.set_proc_list(proc_list);

		for(ImmersedBoundaryDynamics3D<double, DESCRIPTOR, Periodicity>::ObjMapIterator iter = fsi2.particles_begin();
						iter != fsi2.particles_end(); ++iter) {
			plb::fsi::RigidParticle3D<T> * p = dynamic_cast<plb::fsi::RigidParticle3D<T> *>(iter->second);
			for(plint i = 0; i < proc_list.size(); ++i) {
				comm_buff.pack<plint>(proc_list[i], p->get_id());
				comm_buff.pack(proc_list[i], p->center_of_mass());
				comm_buff.pack(proc_list[i], p->orientation());
			}

			positions[p->get_id()] = p->center_of_mass();
			orientations[p->get_id()] = p->orientation();
		}

		comm_buff.send_and_receive_no_wait(true);
		comm_buff.finalize_send_and_receive();

		char * it = comm_buff.recv_buffer_begin();
		while(it != comm_buff.recv_buffer_end()) {
			plint id;
			utils::unpack(it, id);
			utils::unpack(it, positions[id]);
			utils::unpack(it, orientations[id]);
		}

		// Add particles to the immersed boundary object
		{
			int N;
			for(int i = 0, N = 0; i < particleTypes.size(); ++i) {
				for(int j = 0; j < particleTypes[i].numParticles; ++j, ++N) {
					ParticleBase3D<T> * p = particleTypes[i].particle->clone();
					DeformableParticle3D<T> * pd = dynamic_cast<DeformableParticle3D<T> *>(p);
					RigidParticle3D<T> * pr = dynamic_cast<RigidParticle3D<T> *>(p);
					if(pd) {
						Transform<T> transform;
						transform.translate(-pd->center_of_mass())
								 .rotate(orientations[N])
								 .translate(positions[N]);
						pd->transform_vertices(transform);
					} else if(pr) {
						pr->orientation() = orientations[N];
						pr->set_center_of_mass(positions[N]);
					} else {
						std::cerr << "An error occured during particle type casting" << std::endl;
					}
					p->update();

					fsi.add_particle(p);
					delete p;
				}
			}
		}
	}
}



#endif /* CREATE_PARTICLES_H_ */
