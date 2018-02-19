#ifndef BOUNDARY_H_
#define BOUNDARY_H_
#include "geometry.h"
#include <vector>

namespace plb {
namespace fsi {

template<class T>
struct BoundaryNode {
	Array<T, 3> centroid;
	Array<T, 3> normal;
	T area;
};

template<class T>
class Boundary {
public:
	virtual ~Boundary() { }
	virtual bool contains(const Array<T, 3> &, T margin = 0) const = 0;
	virtual bool does_intersect(const geo::Rect<T> &, T margin = 0) const = 0;
	virtual bool distance_to_boundary_less_than(const Array<T, 3> &, T) const = 0;
	virtual T distance_to_boundary(const Array<T, 3> &) const = 0;
	virtual T distance_to_boundary(const Array<T, 3> &, const Array<T, 3> &) const = 0;
	virtual Array<T, 3> get_normal(const Array<T, 3> &) const = 0;
	virtual Box3D get_bounding_box(plint margin) const = 0;
	void get_boundary_nodes(std::vector<BoundaryNode<T> > & res)
	{
		std::vector<Array<T, 3> > vertices;
		std::vector<Array<plint, 3> > triangles;
		get_triangulation(vertices, triangles);

		res.resize(triangles.size());

		for(plint i=0; i < triangles.size(); ++i) {
			Array<T, 3> v0 = vertices[triangles[i][0]];
			Array<T, 3> v1 = vertices[triangles[i][1]];
			Array<T, 3> v2 = vertices[triangles[i][2]];

			// Compute centroid
			res[i].centroid = (v0 + v1 + v2) / (T) 3.0;

			// Get cross product of edges and normal vector.
			Array<T, 3> V1mV0 = v1 - v0;
			Array<T, 3> V2mV0 = v2 - v0;
			Array<T, 3> N = crossProduct(V1mV0, V2mV0);

			// Store normal
			T normal_norm = norm(N);
			res[i].normal = N / normal_norm;

			// Compute area
			res[i].area = 0.5 * normal_norm;
		}
	}

	virtual bool trace_ray(const Array<T, 3> & x0, const Array<T, 3> & x1, T & t) const = 0;

	void writeVTK(std::string filename) {
		if(global::mpi().isMainProcessor()) {
			std::vector<Array<T, 3> > vertices;
			std::vector<Array<plint, 3> > triangles;
			get_triangulation(vertices, triangles);

			std::ofstream out( filename.c_str(), std::ios_base::trunc);

			out << "<?xml version=\"1.0\"?>\n";
			out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
			out << "  <UnstructuredGrid>\n";
			out << "    <Piece NumberOfPoints=\""<< vertices.size() <<"\" NumberOfCells=\""<< triangles.size() <<"\">\n";
			out << "    <Points>\n";
			out << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
			for(plint i = 0; i < vertices.size(); ++i) {
				out << "        " << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << std::endl;
			}
			out << "        </DataArray>\n";
			out << "    </Points>\n";
			out << "    <Cells>\n";
			out << "        <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n";

			for(plint i = 0; i < triangles.size(); ++i) {
				out << "        " << triangles[i][0] << " "
							<< triangles[i][1] << " "
							<< triangles[i][2] << std::endl;
			}
			out << "        </DataArray>\n";
			out << "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n";
			for(plint i=0; i < triangles.size(); i++)
				out << (i+1)*3 << " ";
			out << std::endl;
			out << "        </DataArray>\n";

			out << "        <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
			for(int i=0; i < triangles.size(); i++)
				out << " 5";
			out << std::endl;

			out << "        </DataArray>\n";
			out << "      </Cells>\n";
			out << "    </Piece>\n";
			out << "  </UnstructuredGrid>\n";
			out << "</VTKFile>\n";
			out.close();
		}
	}
private:
	virtual void get_triangulation(std::vector<Array<T, 3> > &, std::vector<Array<plint, 3> > &) const = 0;
};

template<class T>
class BoundaryOutside : public plb::DomainFunctional3D {
public:
	BoundaryOutside(const Boundary<T> & boundary_, plint margin) : boundary(boundary_), margin_(margin) { }

	virtual BoundaryOutside<T> * clone() const { return new BoundaryOutside<T>(*this); }

	virtual bool operator()(plint iX, plint iY, plint iZ) const
	{
		return ! boundary.contains(Array<T, 3>(iX, iY, iZ), (T)-margin_);
	}

private:
	const Boundary<T> & boundary;
	T margin_;
};


/*
 * Pipe boundary aligned with the x-direction
 */
template<class T>
class PipeBoundary : public Boundary<T> {
public:
	PipeBoundary() : x0(0), y0(0), z0(0), radius(0), length(0) { }
	PipeBoundary(T y0, T z0, T radius, T length)
	: x0(0), y0(y0), z0(z0), radius(radius), length(length)
	{ }

	virtual ~PipeBoundary() { }

	Box3D get_bounding_box(plint margin) const
	{
		return Box3D(std::floor(x0), std::ceil(x0+std::floor(length)),
					std::ceil(y0-radius)-margin, std::floor(y0+radius)+margin,
					std::ceil(z0-radius)-margin, std::floor(z0+radius)+margin);
	}

	virtual bool contains(const Array<T,3> & pos, T margin = 0) const
	{
		return (util::sqr(pos[1]-y0) + util::sqr(pos[2]-z0)) < util::sqr(radius - margin);
	}

	virtual bool distance_to_boundary_less_than(const Array<T, 3> & pos, T dist) const
	{
		return (util::sqr(pos[1]-y0) + util::sqr(pos[2]-z0)) > util::sqr(radius-dist);
	}

	virtual T distance_to_boundary(const Array<T, 3> & pos) const
	{
		return radius - std::sqrt(util::sqr(pos[1]-y0) + util::sqr(pos[2]-z0));
	}

	virtual bool trace_ray(const Array<T, 3> & x0, const Array<T, 3> & x1, T & t) const
	{
		const T y1 = x0[1]-y0, y2 = x1[1]-y0, z1 = x0[2]-z0, z2 = x1[2]-z0;

		const T a = (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1);
		const T b = y1*(y2-y1) + z1*(z2-z1);
		const T c = y1*y1 + z1*z1 - radius*radius;

		// Check if a real solution exists
		const T disc = b*b - c*a;
		if(disc < 0. || std::abs(a) < 10*std::numeric_limits<T>::epsilon())
			return false;
		const T bOverA = -b/a;
		const T sqrtDisc = std::sqrt(disc)/a;

		// Test both roots
		t = bOverA - sqrtDisc;
		if(t >= 0. && t <= 1.)
			return true;
		t = bOverA + sqrtDisc;
		if(t >= 0. && t <= 1.)
			return true;
		return false;
	}

	/**
	 * distance_to_boundary
	 *   Params: pos - position to trace from
	 *           dir - direction in which the distance to the boundary is traced
	 *   Returns the smallest distance (in absolute sense) to the boundary
	 */
	virtual T distance_to_boundary(const Array<T, 3> & pos, const Array<T, 3> & dir) const
	{
		// Normalize the direction
		Array<T, 3> d2 = dir / norm(dir);

		const T dy = pos[1] - y0;
		const T dz = pos[2] - z0;
		const T nnrm_sqr = util::sqr(d2[1]) + util::sqr(d2[2]);

		const T a = (d2[1]*dy+ d2[2]*dz) / nnrm_sqr;
		const T b = a*a - (dy*dy + dz*dz - radius*radius)/nnrm_sqr;

		// No solution
		if(b <= 0.0001)
			return std::numeric_limits<T>::max();

		// Find the smallest solution
		const T c = std::sqrt(b);

		const T dist1 = -a + c;
		const T dist2 = -a - c;

		return std::abs(dist1) < std::abs(dist2) ? dist1 : dist2;
	}

	/*virtual void get_mirror_point(const Array<T, 3> & p, Array<T, 3> & res) const
	{
		const T factor = (T) 2 * (radius / std::sqrt(util::sqr(p[1]-y0) + util::sqr(p[2]-z0)) - (T)1);
		res[0] = p[0];
		res[1] = p[1] + factor * (p[1] - y0);
		res[2] = p[2] + factor * (p[2] - z0);
	}*/

	virtual bool does_intersect(const geo::Rect<T> & r, T width = 0) const
	{
		// Calculate the distance between the cylinder center and the closest point in the rectangle
		T tmp = (util::sqr(y0 - geo::clamp(y0, r.y0, r.y1)) +
				util::sqr(z0 - geo::clamp(z0, r.z0, r.z1)));
		// Test if the rectangle intersects the cylinder
		if(tmp <= util::sqr(radius)) {
			// Check if the rectangle is fully contained in the cylinder with radius = radius - width
			// In that case, the rectangle does not intersect the boundary
			return ! (contains(Array<T, 3>((T)0, r.y0, r.z0), width) &&
					  contains(Array<T, 3>((T)0, r.y0, r.z1), width) &&
					  contains(Array<T, 3>((T)0, r.y1, r.z0), width) &&
					  contains(Array<T, 3>((T)0, r.y1, r.z1), width));
		} else {
			return false;
		}
	}

	virtual Array<T, 3> get_normal(const Array<T, 3> & pos) const
	{
		Array<T, 3> d(0, pos[1]-y0, pos[2]-z0);
		return d / norm(d);
	}

	T x0, y0, z0, radius, length;

private:
	virtual void get_triangulation(std::vector<Array<T, 3> > & vertices, std::vector<Array<plint, 3> > & indices) const
	{
		// The elements are scaled so that the resulting triangles have area ~1
		// Number of angular elements
		T perimeter = 2 * M_PI * radius;
		plint ntheta = std::ceil(perimeter / std::sqrt(2));
		T dtheta = 2 * M_PI / ntheta;

		// Number of elements in lengthwise direction
		plint nx = std::ceil((length+1)/std::sqrt(2));
		T dx = length / (nx-1);

		// Number of vertices
		plint nvert = nx*ntheta;

		//std::cout << " Triangulation #vertices = " << nvert << std::endl;

		// Create vertices
		vertices.reserve(nvert);

		for(plint i=0; i < nx; ++i) {
			T theta_start = i*dtheta/2;
			T x = x0 + i*dx;
			for(plint j=0; j < ntheta; ++j) {
				T theta = theta_start + j*dtheta;
				vertices.push_back(Array<T, 3>(x, y0 + radius*std::cos(theta), z0 + radius*std::sin(theta)));
				//std::cout << x << "," << y0 + radius*std::cos(theta) << "," << z0 + radius*std::sin(theta) << std::endl;
			}
		}

		//std::cout << "-----------" << std::endl;

		// Create triangles
		plint n_tri = (nx-1)*2*ntheta;
		indices.reserve(n_tri);
		for(plint i=0; i < (nx-1); ++i) {
			plint ind_base = i*ntheta;
			for(plint j=0; j < (ntheta-1); ++j) {
				indices.push_back(Array<plint, 3>(ind_base+j, ind_base+j+1, ind_base+j+ntheta));
				indices.push_back(Array<plint, 3>(ind_base+j+ntheta, ind_base+j+ntheta+1, ind_base+j+1));
		    }
			indices.push_back(Array<plint, 3>(ind_base+ntheta-1, ind_base+2*ntheta-1, ind_base));
			indices.push_back(Array<plint, 3>(ind_base+2*ntheta-1, ind_base+ntheta, ind_base));
		}
	}
};

/*
 * Couette flow with walls in the xz-plane
 */
template<class T>
class ParallelPlatesBoundary : public Boundary<T> {
public:
	ParallelPlatesBoundary(T y0_, T y1_) : y0(y0_), y1(y1_) { }

	virtual ~ParallelPlatesBoundary() { }

	virtual bool contains(const Array<T, 3> & pos, T margin = 0) const
	{
		return pos[1] > (y0+margin) && pos[1] < (y1-margin);
	}

	virtual bool does_intersect(const geo::Rect<T> & b, T margin = 0) const
	{
		return (b.y0 < (y0 + margin) && b.y1 > (y0 + margin)) ||
				(b.y0 < (y1 - margin) && b.y1 > (y1 - margin));
	}

	virtual bool distance_to_boundary_less_than(const Array<T, 3> & r, T dist) const
	{
		return (r[1] + dist) > y1 || (r[1] - dist) < y0;
	}

	virtual T distance_to_boundary(const Array<T, 3> & r) const
	{
		return (r[1] < 0.5*(y0+y1)) ? r[1] - y0 : y1 - r[1];
	}

	virtual T distance_to_boundary(const Array<T, 3> & r, const Array<T, 3> & dir) const
	{
		// Avoid division by zero
		if(std::abs(dir[1]) < 1e-6)
			return std::numeric_limits<T>::max();

		// Evaluate distance to the boundaries (elementary trigonometry)
		return (r[1] < 0.5*(y0+y1)) ?
				(y0 - r[1]) * norm(dir) / dir[1] : (y1 - r[1]) * norm(dir) / dir[1];
	}

	virtual bool trace_ray(const Array<T, 3> & x0, const Array<T, 3> & x1, T & t) const
	{
		// Upper plate
		t = (y1 - x0[1]) / (x1[1] - x0[1]);
		if(t > 0 && t < 1)
			return true;
		// Lower plate
		t = (y0 - x0[1]) / (x1[1] - x0[1]);
		if(t > 0 && t < 1)
			return true;
		return false;
	}

	virtual Array<T, 3> get_normal(const Array<T, 3> & r) const
	{
		return (r[1] < 0.5*(y0+y1)) ? Array<T, 3>(0, -1.0, 0) : Array<T, 3>(0, 1.0, 0);
	}

	virtual Box3D get_bounding_box(plint margin) const
	{
		return Box3D(0, 0, y0, y1, 0, 0);
	}

private:
	T y0, y1;
	virtual void get_triangulation(std::vector<Array<T, 3> > &, std::vector<Array<plint, 3> > &) const
	{

	}
};

}

}




#endif /* BOUNDARY_H_ */
