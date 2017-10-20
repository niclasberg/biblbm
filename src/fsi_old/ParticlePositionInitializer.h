#ifndef PARTICLEPOSITIONINITIALIZER_H_
#define PARTICLEPOSITIONINITIALIZER_H_
#include "Grid.h"
#include <vector>
#include "core/globalDefs.h"
#include "geometry.h"

namespace plb {

namespace fsi {

template<class T>
struct ParticlePosition {
	char tag;
	Array<T, 3> pos;
	T radius;
};

template<class T>
class ParticlePositionInitializer {
public:
	ParticlePositionInitializer(const Box3D & bounding_box_, Boundary<T> * boundary_, T max_radius)
	: grid(geo::Rect<T>(bounding_box_), 2*max_radius),  bounding_box(bounding_box_), boundary(boundary_)
	{

	}

	bool generate_points(pluint N, T radius, char tag);
	void get_points(char tag, std::vector<Array<T, 3> > & points);

private:
	bool generate_position(ParticlePosition<T> &);

	std::vector<ParticlePosition<T> > positions;
	Grid<T, ParticlePosition<T> > grid;
	Boundary<T> * boundary;
	Box3D bounding_box;
};

} /* namespace fsi */

} /* namespace plb */


#endif /* PARTICLEPOSITIONINITIALIZER_H_ */
