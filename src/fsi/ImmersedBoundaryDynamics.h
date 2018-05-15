/*
 * ImmersedBoundaryDynamics.h
 *
 *  Created on: 22 feb 2014
 *      Author: niber
 */

#ifndef IMMERSEDBOUNDARYDYNAMICS_H_
#define IMMERSEDBOUNDARYDYNAMICS_H_
#include "ParticleBase.h"
#include "Dirac.h"
#include "CommunicationBuffer.h"
#include <vector>
#include <map>
#include "Grid.h"
#include "Periodicity.h"
#include "geometry.h"
#include "ZBuffer.h"
#include "ForceDecorator.h"
#include "Boundary.h"
#include "ParticleParticleInteraction.h"

namespace plb {

namespace fsi {

// Forward declarations
template<class T> class Boundary;
template<class T> class ParticleShapeLibrary;

template<class T>
struct NonLocalNode {
	plint proc_id;
	plint node_id;
	plint particle_id;
	Array<T, 3> pos;
	Array<T, 3> vel;
	Array<T, 3> force;
};

template<class T>
struct MpiCommInfo {
	plint proc_id;
	plint num_nodes_to_receive;
	geo::Rect<T> domain_with_fsi_envelope;
	Box3D domain;
};

/* Main class for the Immersed boundary interaction */
template<class T, template<typename U> class Descriptor, class Periodicity>
class ImmersedBoundaryDynamics3D {
private:
	// Typedefs
	//typedef TopHatDirac<T, 3> Dirac;
	typedef RomaDirac<T, 3> Dirac;

public:
	typedef ParticleBase3D<T> ParticleType;
	typedef typename std::map<plint, ParticleType * > ObjMapType;
	typedef typename ObjMapType::iterator ObjMapIterator;
	typedef typename ObjMapType::const_iterator ObjMapConstIterator;

	ImmersedBoundaryDynamics3D(const MultiBlockLattice3D<T, Descriptor> &, const ParticleShapeLibrary<T> &);
	~ImmersedBoundaryDynamics3D();

	// Initialization (should be called after all particles have been added)
	void init(const Box3D &, BlockLattice3D<T, Descriptor> &);
	void init();
	bool is_init() const { return is_init_; }

	// Particle state synchronization
	void synchronize_particle_states();
	void synchronize_particle_states_and_voxelization(const Box3D &, BlockLattice3D<T, Descriptor> &);

	void clear();

	// Force decorators
	void add_force_decorator(ForceDecorator<T> * force_decorator) { force_decorators_.push_back(force_decorator); }
	void add_pp_force(PPInteractionForce<T> * pp_force) { pp_forces_.push_back(pp_force); }

	// Fsi methods
	void set_forces_to_zero();
	void interpolate_velocity(const Box3D &, TensorField3D<T, 3> &);
	void compute_and_spread_forces(const Box3D &, BlockLattice3D<T, Descriptor> &);
	void move_vertices_and_revoxelize(const Box3D &, BlockLattice3D<T, Descriptor> &);
	void move_vertices();
	template<class InteractionFunctional> void compute_collision_forces(const InteractionFunctional &, const Boundary<T> * = 0);

	// Add a particle to the simulation. Only particles contained in the domain of
	// this processor are saved.
	void add_particle(const ParticleType *);

	void set_boundary(Boundary<T> * boundary_) { this->boundary = boundary_; }

	// Particle access
	ObjMapIterator particles_begin() { return particles.begin(); }
	ObjMapIterator particles_end() { return particles.end(); }
	pluint count_particles() const { return particles.size(); }
	bool get_particle(plint i, ParticleType *& particle);
	bool get_particle(plint i, const ParticleType *& particle) const;

	// Global reductions
	plint count_total_particles() const;
	T compute_total_volume() const;
	void get_all_particles(std::vector<const ParticleType *> &) const;

	// Output
	void write_particles_as_vtk(pluint it);
	void write_lightweight_particle_data(pluint it);

	// Checkpoint files
	void save_checkpoint(FileName) const;
	void load_checkpoint(FileName);

private:
	void add_particle(const ParticleType *, plint id);
	void add_nonlocal_node(NonLocalNode<T> &);

	// Mpi communication initialization
	void initialize_domains();
	void find_blocks_in_domain(const geo::Rect<T> & bb, std::vector<plint> & neighboring_blocks) const;

	// Velocity interpolation methods
	void send_interpolation_data();
	void receive_and_reduce_interpolation_data();

	// Force computation and spreading
	void send_force_data();
	void receive_force_data();
	void compute_forces();
	void compute_collision_forces();
	template<class InteractionFunctional> void compute_particle_particle_interaction_impl(const InteractionFunctional &, const Boundary<T> *);
	void pack_local_particle_forces();
	void unpack_local_particle_forces(char *&);
	void pack_nonlocal_node_forces();
	void reduce_nonlocal_node_forces(char *&);

	// Vertex movement
	void move_vertices_impl();
	void voxelize();
	void apply_voxelization(const Box3D &, BlockLattice3D<T, Descriptor> &);
	void create_nodes();
	void create_nodes_from_particle(ParticleType &);
	void clear_nodes();
	void add_nodes_from_particle(const ParticleType &);

	// Mpi communication methods
	void pack_handoff_particles(std::vector<plint> &);
	void unpack_handoff_particles(char *&);
	void pack_local_nodes();
	void unpack_nonlocal_nodes(char *&);
	void pack_voxelization_data();
	void unpack_voxelization_data(char *&);

	// Force spreading
	template<class NodeType>
	void spread_forces_impl(const Box3D &, BlockLattice3D<T, Descriptor>  &, std::vector<NodeType> &);

	//
	bool intersects(const ParticleType &) const;
	bool is_master_of(const ParticleType &) const;
	plint get_proc_id(const ParticleType &) const;

	// Id generation for particles
	pluint get_next_id() { return num_particles++; }

	geo::Rect<T> local_particle_bounding_box() const;

private:
	Boundary<T> * boundary;

	// Mpi communication
	CommunicationBuffer comm_buffer;
	std::map<plint, MpiCommInfo<T> > comm_info;

	// Optimized evaluation of the dirac function
	SampledDirac<T, Dirac, 3> sampled_dirac;

	// Node containers
	std::vector<Vertex<T> *> local_nodes;
	std::vector<Vertex<T> *> local_nodes_envelope;
	std::vector<Vertex<T> *> local_nodes_boundary;
	std::vector<NonLocalNode<T> > nonlocal_nodes;
	std::vector<NonLocalNode<T> > nonlocal_nodes_envelope;
	std::vector<NonLocalNode<T> > nonlocal_nodes_boundary;
	std::map<plint, plint> particle_start_id_;

	// Geometry arithmetic (contains all methods to implement periodicity)
	typedef typename Periodicity::ArithmeticType ArithmeticType;
	ArithmeticType arithmetic;

	// Library of all particle shapes
	const ParticleShapeLibrary<T> & shape_library;

	// Voxelization
	ZBuffer<T, typename Periodicity::ArithmeticType> z_buffer;

	// Domain data
	const MultiBlockManagement3D & management;	// Underlying block management
	Box3D global_bounding_box;					// Global bounding box of the domain
	pluint nx, ny, nz;							// Global domain dimensions
	Box3D domain;								// Domain for the current MPI process
	geo::Rect<T> domain_with_fsi_envelope;		// Domain including the envelope for fsi communication
	geo::Rect<T> bulk_domain;					// Domain excluding the envelope (i.e. the domain with
												// 	fluid nodes with no data dependencies to other processors)

	ObjMapType particles;						// Particles whose bounding sphere intersects with the current domain
	pluint num_particles;						// Number of particles (also the maximal particle id + 1)
	plint num_nodes;
	std::vector<plint> neighbouring_proc_ids;

	// Force decorators
	std::vector<ForceDecorator<T> *> force_decorators_;
	std::vector<PPInteractionForce<T> *> pp_forces_;

	bool is_init_;
};

} /* namespace fsi */

} /* namespace plb */

#endif /* IMMERSEDBOUNDARYDYNAMICS_H_ */
