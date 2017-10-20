#include "palabos3D.h"
#include "palabos3D.hh"
#include "fsi/headers.h"
#include "fsi/headers.hh"
#include "fsi/Potentials.h"
#include "../common/rbc_parameters.h"
#include "../common/calibrate_rbc_shape.h"

using namespace plb;
using namespace fsi;

#define DESCRIPTOR descriptors::ForcedPhaseD3Q19Descriptor
typedef Periodicity3D<double, true, false, false> Periodicity;

int main(int argc, char * argv[])
{
	plbInit(&argc, &argv);

	// Create boundary
	PipeBoundary<double> boundary(10);

	// Create lattice
	MultiBlockLattice3D<double, DESCRIPTOR> lattice;
	Box3D bb = lattice.getBoundingBox();

	// Generate RBC shape
	ParticleShapeLibrary<double> shape_library;
	shape_library.read_and_store_mesh("sphere_TRImesh_626n_1248el.msh", "RBC");
	ParticleShape<double> * rbc_shape = shape_library.get_by_tag("RBC");

	// Create fsi
	ImmersedBoundaryDynamics3D<double, DESCRIPTOR, Periodicity> fsi(lattice, shape_library);

	// Create RBC and calibrate shape
	RBCParticle<double> rbc(rbc_shape, create_rbc_params(rbc_shape, 0.01, 0.1, 0.001));
	shrink_rbc_volume(rbc, 0.59*rbc_shape->get_volume(), 4e4);
	rbc.set_center_of_mass(Array<double, 3>(0, 0, 0));
	rbc.set_minor_axis_orientation(Array<double, 3>(0, 0, 1));

	// Convert to shape representation
	rbc.store_shape(shape_library, "SOLID_RBC");
	ParticleShape<double> * solid_rbc_shape = shape_library.get_by_tag("SOLID_RBC");

	// Create particle
	plb::fsi::RigidParticle3D<double> rbc_rigid(solid_rbc_shape);
	rbc_rigid.scale() = 0.2;
	rbc_rigid.update();

	// Generate initial distribution
	std::vector<Array<double, 3> > positions;
	plint N_rbc = 100;
	if(global::mpi().isMainProcessor()) {
		double rbc_radius = rbc_rigid.get_radius() + 1.;
		ParticlePositionInitializer<double> position_initializer(bb, &boundary, rbc_radius);
		position_initializer.generate_points(N_rbc, rbc_radius, 'R');
		position_initializer.get_points('R', positions);
	} else {
		positions.resize(N_rbc);
	}

	// Broadcast
	global::mpi().bCast((double *)positions.data(), N_rbc*3);

	// Create and insert particles to the fsi object
	for(plint i = 0; i < N_rbc; ++i) {
		rbc_rigid.set_center_of_mass(positions[i]);
		fsi.add_particle(&rbc_rigid);
	}

	// Initialize fsi
	fsi.init();

	// Potential force
	SpringPotential<double> potential(1., 1e-5);

	// Grow the particles while under the influence of collision forces
	double scale;
	for(plint it = 0; it < Nit; ++it) {
		if((it % 10) == 0) {
			// Rescale all particles
			double fraction = (double)it / (double) (Nit - 1);
			scale = initial_scale + (1. - initial_scale)*fraction;

			for(ImmersedBoundaryDynamics3D<double, DESCRIPTOR, Periodicity>::ObjMapIterator iter = fsi.particles_begin();
					iter != fsi.particles_end(); ++iter) {
				RigidParticle3D<double> * p = dynamic_cast<RigidParticle3D<double> *>(iter->second);
				if(p) {
					p->scale() = scale;
					p->update();
				}
			}

			fsi.synchronize_particle_states();
		}

		fsi.compute_collision_forces(potential);
		fsi.move_vertices();
	}
}
