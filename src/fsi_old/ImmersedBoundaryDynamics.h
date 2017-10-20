/*
 * ImmersedBoundaryDynamics.h
 *
 *  Created on: 22 feb 2014
 *      Author: niber
 */

#ifndef IMMERSEDBOUNDARYDYNAMICS_H_
#define IMMERSEDBOUNDARYDYNAMICS_H_
#include "RigidParticle.h"
#include "Dirac.h"
#include "CommunicationBuffer.h"
#include <vector>
#include <map>
#include "Grid.h"
#include "Periodicity.h"
#include "Quadrature.h"

namespace plb {

namespace fsi {

// Forward declarations
template<class T> class Boundary;
template<class T> class ParticleShapeLibrary;

template<class T>
struct SolidNode {
	Array<T, 3> vel; 		// Solid body velocity at the node
	Array<T, 3> pos;		// Node position (world frame)
	Array<T, 3> pos_rel;	// Node position (relative to the body center)
	Array<T, 3> normal;		// Element normal
	T volume;					// Area of the surface element
	plint obj_id, node_id;	// Object and node ids
	RigidParticle3D<T> * particle;
};

// Solid node container
// This class will never deallocate memory when cleared
// and memory allocations will thus be kept at minimum.
template<class T>
class SolidNodeList {
public:
	SolidNodeList() : solid_node_capacity(1), solid_node_count(0), solid_nodes(new SolidNode<T>[1]) { }
	virtual ~SolidNodeList() { delete [] solid_nodes; }
	SolidNode<T> & operator[](pluint i) { PLB_PRECONDITION(i < solid_node_count); return solid_nodes[i]; }
	const SolidNode<T> & operator[](pluint i) const { PLB_PRECONDITION(i < solid_node_count); return solid_nodes[i]; }
	const pluint & size() const { return solid_node_count; }
	void clear() { solid_node_count = 0; }
	SolidNode<T> & last() { return solid_nodes[solid_node_count-1]; }
	void push_node(const SolidNode<T> & node);

private:
	// Disable copy constructor and assignment operator
	SolidNodeList(const SolidNodeList &);
	SolidNodeList & operator=(const SolidNodeList &);

	SolidNode<T> * solid_nodes;
	pluint solid_node_capacity, solid_node_count;
};

/* Main class for the Immersed boundary interaction */
template<class T, template<typename U> class Descriptor, class Periodicity>
class ImmersedBoundaryDynamics3D {
private:
	// Typedefs
	typedef RomaDirac<T, 3> Dirac;
	typedef Grid<T, SolidNode<T> > GridType;
	typedef std::vector<typename GridType::CellType *> BoundaryCellContainerType;
	typedef FsiForceCommunicator<T, Dirac, typename Periodicity::ArithmeticType> FsiForceCommunicatorType;
	typedef ParticleForceCommunicator<T, typename Periodicity::ArithmeticType> ParticleForceCommunicatorType;
	typedef ParticleStateCommunicator<T, typename Periodicity::ArithmeticType> ParticleStateCommunicatorType;

public:
	typedef typename std::map<plint, RigidParticle3D<T> * > ObjMapType;
	typedef typename ObjMapType::iterator ObjMapIterator;

	ImmersedBoundaryDynamics3D(const MultiBlockLattice3D<T, Descriptor> &, const ParticleShapeLibrary<T> &);
	~ImmersedBoundaryDynamics3D();

	/* This method does the following
	 * 1. Moves all particles based on the current linear and angular velocities
	 * 2. Computes particle-particle and wall-particle interaction forces
	 * 3. Computes the fluid-structure interaction forces and updates the linear and angular momentum of the particles
	 * The 3rd step can be set to be repeated multiple times by changing the value of fsi_subiterations.
	 *
	 * NOTE: This method should not be called explicitly, but rather implicitly through the
	 * ImmersedBoundaryWrapperFunctional3D data processor 														*/
	template<class InteractionFunctional>
	void do_fsi(
			Box3D,
			BlockLattice3D<T, Descriptor> &,
			ScalarField3D<T> &,
			TensorField3D<T, 3> &,
			const Array<T, 3> &,
			const InteractionFunctional &);

	// Dry iteration (no fsi, only integration of particle equations of motion under the effect of collision forces)
	// This is useful for creating initial conditions.
	template<class InteractionFunctional>
	void do_dry_iteration(const InteractionFunctional &, T damping = 0.3);

	// Add a particle to the simulation. Only particles contained in the domain of
	// this processor are saved.
	void add_particle(RigidParticle3D<T> *);

	// Total number of particles (spread over all mpi processes)
	// NOTE: this method performs a global reduction and is slow. Use with care.
	pluint count_particles() const;

	// Check if we have overlapping particles
	// NOTE: This method is SLOW and should not be used in any critical loops
	// TODO: implement
	bool has_overlapping_particles() const;

	// Access particle iterators
	ObjMapIterator particles_begin() { return particles.begin(); }
	ObjMapIterator particles_end() { return particles.end(); }

	// Setters and getters for the number of fsi subiterations
	void set_fsi_subiteration_count(pluint value) { fsi_subiterations = value; }
	pluint get_fsi_subiteration_count() const { return fsi_subiterations; }

	// Collision interaction distance (no collision forces are not computed if the
	// distance between two nodes is greater than this value).
	// NOTE: this must be set before particles are added to the domain
	void set_interaction_cutoff_distance(T);
	T get_interaction_cutoff_distance() const { return interaction_distance; }

	// Set the shape of the boundary. Can be set to NULL to remove the boundary
	void set_boundary(Boundary<T> * boundary) { this->boundary = boundary; }
	bool has_boundary() const { return boundary; }

	// Output
	void write_particles_as_vtk(pluint it);
	void write_particle_states(pluint it);

	// set current iteration
	void set_iteration(plint iteration_) { iteration = iteration_; }

	// Checkpoint files
	void save_checkpoint(FileName) const;
	void load_checkpoint(FileName);

private:
	// Move particles and communicate new states with neighboring processors
	void move_particles(T damping = 0);

	// Interaction force evaluation
	template<class InteractionFunctional>
	void compute_interaction_forces(const InteractionFunctional &);
	template<class InteractionFunctional>
	void collide_all_in_cells(typename GridType::CellType & cell1, typename GridType::CellType & cell2, const InteractionFunctional &);
	template<class InteractionFunctional>
	void handle_collision(const SolidNode<T> &, const SolidNode<T> &, const InteractionFunctional &);

	void test_collision_and_append_on_success(const SolidNode<T> &, SolidNode<T> &, std::vector<SolidNode<T> *> &);

	// Compute fsi forces on fluid and particles
	void compute_fsi_forces(Box3D, BlockLattice3D<T, Descriptor> &, ScalarField3D<T> &, TensorField3D<T, 3> &);
	void compute_envelope_fsi_forces(Box3D, BlockLattice3D<T, Descriptor> &, ScalarField3D<T> &, TensorField3D<T, 3> &);
	void compute_bulk_fsi_forces(Box3D, BlockLattice3D<T, Descriptor> &, ScalarField3D<T> &, TensorField3D<T, 3> &);
	void spread_envelope_fsi_forces(Box3D, BlockLattice3D<T, Descriptor> &, ScalarField3D<T> &, TensorField3D<T, 3> &);
	void update_fluid_velocity(Box3D, BlockLattice3D<T, Descriptor> &, TensorField3D<T, 3> &);

	// Loops over the mesh of a particle and fills the local_node_list, envelope_node_list and the grid
	// with SolidNodes.
	void add_solid_nodes_from_particle(RigidParticle3D<T> &);

	// Recompute the velocity of each node
	void update_solid_nodes_velocities();

	// Id generation for particles
	pluint get_next_id() { return num_particles++; }

	void set_particle_proc_id(RigidParticle3D<T> &);

	// These methods are responsible for setting up the communicators with processor-domain information
	void rebuild_neighbor_proc_list();
	void find_blocks_in_domain(const geo::Rect<T> &, std::vector<plint> &) const;
	void initialize_fsi_communicator();
	void initialize_periodicity(const MultiBlockLattice3D<T, Descriptor> &);
	void initialize_domains();
	void create_proc_map(const geo::Rect<T> &, T, std::vector<plint> &, std::vector<std::pair<geo::Rect<T>, plint> > &) const;

private:
	// Optimized evaluation of the dirac function
	SampledDirac<T, Dirac, 3> sampled_dirac;

	// Geometry arithmetic (contains all methods to implement periodicity)
	typename Periodicity::ArithmeticType arithmetic;

	// Library of all particle shapes
	const ParticleShapeLibrary<T> & shape_library;

	// Domain data
	const MultiBlockManagement3D & management;	// Underlying block management
	Box3D global_bounding_box;					// Global bounding box of the domain
	pluint nx, ny, nz;							// Global domain dimensions
	Box3D domain;								// Domain for the current MPI process
	geo::Rect<T> domain_with_fsi_envelope;		// Domain including the envelope for fsi communication
	geo::Rect<T> bulk_domain;					// Domain excluding the envelope (i.e. the domain with
												// 	fluid nodes with no data dependencies to other processors)
	geo::Rect<T> interaction_domain;			// Domain including the envelope for interaction computations
	geo::Rect<T> domain_bounds;					// Domain bounding box including the interaction and fsi envelopes
	bool has_periodic_edge[3];					// Flags to set if the domain has edges along which the grid is periodic

	// Node lists and object maps
	SolidNodeList<T> local_node_list;			// Nodes contained in the bulk_domain
	SolidNodeList<T> envelope_node_list;		// Nodes in the envelope. These nodes require special treatment since they
												// 	depend on entities in other MPI processes.

	T interaction_distance;						// Cutoff distance for particle-particle interactions
	T interaction_distance_sqr;					// Cutoff distance squared
	ObjMapType particles;						// Particles whose bounding sphere intersects with the current domain
	T particle_max_radius;						// Max radius of all particles
	pluint num_particles;						// Number of particles (also the maximal particle id + 1)

	pluint fsi_subiterations;					// Number of fsi subiterations
	T acceleration_factor;						// Parameter controlling the convergence of the fsi iterations

	// Communicators
	FsiForceCommunicatorType fsi_force_communicator; 			// Communicator for velocity interpolation
	ParticleForceCommunicatorType particle_force_communicator;	// Communicator for particle force reductions
	ParticleStateCommunicatorType particle_state_communicator;	// Communicator for particle states (momentum, position etc)

	Boundary<T> * boundary;						// Boundary shape. Only used to evaluate particle-wall interactions
	BoundaryCellContainerType boundary_cells;	// List of all Grid cells adjacent to the boundary

	plint iteration;
	bool is_init;

	std::ostream * output_stream;				// Output stream for lightweight particle data
};

template<class T, template<typename U> class Descriptor, class Periodicity, class InteractionFunctional>
class ImmersedBoundaryWrapperFunctional3D : public BoxProcessingFunctional3D {
public:
	ImmersedBoundaryWrapperFunctional3D(
			ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity> & ibm_,
			const Array<T, 3> & force,
			const InteractionFunctional & interaction_)
	: ibm(ibm_),
	  external_force(force),
	  interaction(interaction_)
	{ }

	virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks)
	{
		ibm.do_fsi(domain,
				*(dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[0])),
				*(dynamic_cast<ScalarField3D<T> *>(blocks[1])),
				*(dynamic_cast<TensorField3D<T, 3> *>(blocks[2])),
				external_force,
				interaction
				);
	}

	virtual ImmersedBoundaryWrapperFunctional3D * clone() const
	{
		return new ImmersedBoundaryWrapperFunctional3D(*this);
	}

	virtual BlockDomain::DomainT appliesTo() const {
		return BlockDomain::bulk;
	}

	virtual void getModificationPattern(std::vector<bool>& isWritten) const {
		isWritten[0] = true;
		isWritten[1] = false;
		isWritten[2] = false;
	}

	virtual void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
		modified[0] = modif::staticVariables;
		modified[1] = modif::nothing;
		modified[2] = modif::nothing;
	}
private:
	ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity> & ibm;
	Array<T, 3> external_force;
	InteractionFunctional interaction;
};

template<class T, template<typename U> class Descriptor, class Periodicity, class InteractionFunctional>
ImmersedBoundaryWrapperFunctional3D<T, Descriptor, Periodicity, InteractionFunctional> *
wrap_ibm_dynamics3D(ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity> & ibm, const InteractionFunctional & interaction)
{
	return new ImmersedBoundaryWrapperFunctional3D<T, Descriptor, Periodicity, InteractionFunctional>(
			ibm,
			Array<T, 3>(0, 0, 0),
			interaction
		);
}

template<class T, template<typename U> class Descriptor, class Periodicity, class InteractionFunctional>
ImmersedBoundaryWrapperFunctional3D<T, Descriptor, Periodicity, InteractionFunctional> *
wrap_ibm_dynamics3D(ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity> & ibm, const Array<T, 3> & force, const InteractionFunctional & interaction)
{
	return new ImmersedBoundaryWrapperFunctional3D<T, Descriptor, Periodicity, InteractionFunctional>(
			ibm,
			force,
			interaction
		);
}

} /* namespace fsi */

} /* namespace plb */

#endif /* IMMERSEDBOUNDARYDYNAMICS_H_ */
