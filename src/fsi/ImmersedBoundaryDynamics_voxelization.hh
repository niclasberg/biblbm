/*
 * ImmersedBoundaryDynamics_voxelization.hh
 *
 *  Created on: Apr 20, 2016
 *      Author: niber
 */

#ifndef IMMERSEDBOUNDARYDYNAMICS_VOXELIZATION_HH_
#define IMMERSEDBOUNDARYDYNAMICS_VOXELIZATION_HH_
#include "ImmersedBoundaryDynamics.h"

namespace plb {

namespace fsi {

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::voxelize()
{
	// Set z_buffer bounding box
	geo::Rect<T> bb = local_particle_bounding_box();
	z_buffer.set_extent(std::floor(bb.x0), std::ceil(bb.x1), std::floor(bb.y0), std::ceil(bb.y1));

	// Clear data
	z_buffer.clear();

	for(ObjMapIterator it = particles.begin(); it != particles.end(); ++it)
		it->second->voxelizer().find_intersections_z(z_buffer);
}

template<class T, template<typename U> class Descriptor, class Periodicity>
geo::Rect<T> ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::local_particle_bounding_box() const
{
	// Create a bounding box that encompasses the bulk domain and all local particles
	geo::Rect<T> bb = domain;
	for(typename ObjMapType::const_iterator it = particles.begin(); it != particles.end(); ++it) {
		geo::Rect<T> bb_p = it->second->bounding_box();
		bb.x0 = std::min(bb.x0, bb_p.x0);
		bb.x1 = std::max(bb.x1, bb_p.x1);
		bb.y0 = std::min(bb.y0, bb_p.y0);
		bb.y1 = std::max(bb.y1, bb_p.y1);
		bb.z0 = std::min(bb.z0, bb_p.z0);
		bb.z1 = std::max(bb.z1, bb_p.z1);
	}

	// Small fix to avoid creating to large bounding boxes
	if((bb.x1 - bb.x0 + 1) > global_bounding_box.getNx()) {
		bb.x0 = global_bounding_box.x0;
		bb.x1 = global_bounding_box.x1;
	}

	if((bb.y1 - bb.y0 + 1) > global_bounding_box.getNy()) {
		bb.y0 = global_bounding_box.y0;
		bb.y1 = global_bounding_box.y1;
	}

	if((bb.z1 - bb.z0 + 1) > global_bounding_box.getNz()) {
		bb.z0 = global_bounding_box.z0;
		bb.z1 = global_bounding_box.z1;
	}

	return bb;
}

/* Packing of voxelization data */
namespace detail {

template<class T, class Arithmetic>
class ZBufferPacker {
public:
	ZBufferPacker(const ZBuffer<T, Arithmetic> & parent, plint x0, plint x1, plint y0, plint y1)
	: parent_(parent), x0_(x0), x1_(x1), y0_(y0), y1_(y1)
	  { }

	template<class BufferType>
	void pack(BufferType & buf) const {
		parent_.pack(x0_, x1_, y0_, y1_, buf);
	}
private:
	plint x0_, x1_, y0_, y1_;
	const ZBuffer<T, Arithmetic> & parent_;
};

}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::pack_voxelization_data()
{
	for(typename std::map<plint, MpiCommInfo<T> >::iterator it = comm_info.begin(); it != comm_info.end(); ++it) {
		comm_buffer.delegate_pack(it->first,
				detail::ZBufferPacker<T, typename Periodicity::ArithmeticType>(z_buffer,
						it->second.domain.x0, it->second.domain.x1,
						it->second.domain.y0, it->second.domain.y1));
	}
}

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::unpack_voxelization_data(char *& it)
{
	z_buffer.unpack(it);
}

namespace detail {

template<class T, template<typename U> class Descriptor>
inline void set_phase(plint iX, plint iY, plint iZ, BlockLattice3D<T, Descriptor> & lattice, T value)
{
	*(lattice.get(iX, iY, iZ).getExternal(Descriptor<T>::ExternalField::phaseBeginsAt)) = value;
}

template<class T, template<typename U> class Descriptor>
inline void set_phase(plint iX, plint iY, plint iZmin, plint iZmax, BlockLattice3D<T, Descriptor> & lattice, T value)
{
	for(plint iZ = iZmin; iZ <= iZmax; ++iZ)
		set_phase(iX, iY, iZ, lattice, value);
}

} /* namespace detail */

template<class T, template<typename U> class Descriptor, class Periodicity>
void ImmersedBoundaryDynamics3D<T, Descriptor, Periodicity>::apply_voxelization(const Box3D & dom, BlockLattice3D<T, Descriptor> & lattice)
{
	typedef typename ZBuffer<T, typename Periodicity::ArithmeticType>::const_iterator ZBufferIterator;

	Dot3D offset = lattice.getLocation();

	// Sort the z-buffer according to z-value
	z_buffer.sort_entries();

	for(plint iX = dom.x0; iX <= dom.x1; ++iX) {
		for(plint iY = dom.y0; iY <= dom.y1; ++iY) {
			// Get iterators for the range of z-values (smallest to largest) at which
			// a particle surface crosses the z-axis at (iX, iY).
			std::pair<ZBufferIterator, ZBufferIterator> its
				= z_buffer.get(iX+offset.x, iY+offset.y);

			// Fill the regions between each consecutive particle surface intersection
			// The idea is as follows: We start at a point outside the particles so
			// all nodes between this point and the first particle crossing are
			// fluid nodes. Continuing from this point, all nodes until the next crossing
			// are inside nodes. This is process is repeated until the whole domain
			// spanned by this processor is filled.
			bool is_inside = false;
			T z_last = -10000;
			for(ZBufferIterator it = its.first; it != its.second; ++it) {
				detail::set_phase(iX, iY,
						std::max((plint) std::ceil(z_last), domain.z0)-offset.z,
						std::min((plint) std::floor(*it), domain.z1)-offset.z,
						lattice,
						is_inside ? 1.0 : 0.0);
				z_last = *it;
				is_inside = !is_inside;
			}

			// Fill remaining region (above last z-value)
			if(std::ceil(z_last) < domain.z1) {
				detail::set_phase(iX, iY,
						std::max((plint) std::ceil(z_last), domain.z0)-offset.z,
						domain.z1-offset.z,
						lattice,
						is_inside ? 1.0 : 0.0);
			}
		}
	}
}

}

}


#endif /* IMMERSEDBOUNDARYDYNAMICS_VOXELIZATION_HH_ */
