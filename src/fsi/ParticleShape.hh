/*
 * ParticleShape.hh
 *
 *  Created on: 25 dec 2013
 *      Author: niclas
 */

#ifndef PARTICLESHAPE_HH_
#define PARTICLESHAPE_HH_
#include "ParticleShape.h"
#include <limits>
#include "core/util.h"
#include <set>
#include "TriangleUtils.h"

namespace plb {

namespace fsi {

template<class T>
ParticleShape<T>::ParticleShape(Array<T, 3> * verts, unsigned int * indices, unsigned int num_vertices, unsigned int num_elems)
: vertices_(verts, verts+num_vertices), area(0.), volume(0.), center(0, 0, 0)
{
	// Create triangles
	const unsigned int * inds = indices;
	for(plint i = 0; i < num_elems; ++i) {
		Triangle<T> t;
		t.i0 = *(inds++);
		t.i1 = *(inds++);
		t.i2 = *(inds++);

		// Compute properties
		Array<T, 3> normal;
		tri::normal_and_area(vertices_[t.i0], vertices_[t.i1], vertices_[t.i2], normal, t.area);
		t.l_mean = (norm(vertices_[t.i1]-vertices_[t.i0]) + norm(vertices_[t.i2]-vertices_[t.i1]) + norm(vertices_[t.i0] - vertices_[t.i2])) / 3.0;

		area += t.area;

		Array<T, 3> centroid;
		tri::centroid(vertices_[t.i0], vertices_[t.i1], vertices_[t.i2], centroid);

		volume += dot(centroid, normal) * t.area / 3.0;
		triangles_.push_back(t);

		center += centroid * t.area;
	}

	compute_properties();

	center /= area;

	// Find links and adjacent triangles
	compute_connectivity();

	// Compute bounding volumes (box and sphere)
	compute_bounding_volumes();
}

namespace detail {
	template<class T>
	struct LinkCompare {
		bool operator()(const Link<T> & l1, const Link<T> & l2) const
		{
			return l1.i0 == l2.i0 ? l1.i1 < l2.i1 : l1.i0 < l2.i0;
		}
	};

	template<class T>
	void permute_triangle_vertices(Triangle<T> & t) {
		plint tmp = t.i0;
		t.i0 = t.i2;
		t.i2 = t.i1;
		t.i1 = tmp;
	}

	template<class T>
	plint find_missing_vertex(plint i0, plint i1, const Triangle<T> & t)
	{
		if((t.i0 == i0 && t.i1 == i1) || (t.i0 == i1 && t.i1 == i0))
			return t.i2;
		if((t.i1 == i0 && t.i2 == i1) || (t.i1 == i1 && t.i2 == i0))
			return t.i0;
		if((t.i2 == i0 && t.i0 == i1) || (t.i2 == i1 && t.i0 == i0))
			return t.i1;
		return -1;
	}
}

template<class T>
void ParticleShape<T>::compute_connectivity()
{
	// Use a set to store all links for now. This allows for
	// rapid duplicate checking (O(log(N)))
	typedef std::set<Link<T>, detail::LinkCompare<T> > link_set;
	link_set links_found;

	// Go through all triangles and find all unique links
	plint inds[4];
	Link<T> test_link;

	plint num_adj_tris = 0;

	for(plint i = 0; i < triangles_.size(); ++i) {
		// Extract indices to iterate over
		inds[0] = triangles_[i].i0;
		inds[1] = triangles_[i].i1;
		inds[2] = triangles_[i].i2;
		inds[3] = triangles_[i].i0;

		for(plint j = 0; j < 3; ++j) {
			test_link.i0 = std::min(inds[j], inds[j+1]);
			test_link.i1 = std::max(inds[j], inds[j+1]);
			test_link.t1 = -1;
			test_link.t2 = -1;
			test_link.length = norm(vertices_[test_link.i0] - vertices_[test_link.i1]);

			//std::cout << "Link (" << test_link.i0 << ", " << test_link.i1 << "), length: " << test_link.length << std::endl;

			// Check if the link already has been found
			typename link_set::iterator it = links_found.find(test_link);
			if(it == links_found.end()) {
				// Save index of the triangle
				test_link.t1 = i;

				// Save
				links_found.insert(test_link);
			} else {
				// The link already exists, set the second triangle index to the current one
				Link<T> & l2 = const_cast<Link<T> &>(*it);

				if(l2.t1 != i) {
					l2.t2 = i;
					++num_adj_tris;
				}
			}
		}
	}

	// Store the links in an array now, this container is more rapid to iterate over
	links.clear();
	links.insert(links.begin(), links_found.begin(), links_found.end());
	links_found.clear();

	// Create adjacent triangles
	// These are defined by 4 indices, such that the second and third are shared between the
	// triangles and the (0, 1, 2) and (0, 2, 3) vertices makes up the two triangles.
	std::map<plint, std::vector<plint> > neighbors;
	Triangle<T> t1, t2;
	plint i0, i1;
	AdjacentTrianglePair<T> tri_pair;
	for(plint i = 0; i < links.size(); ++i) {
		// Check that the link has two adjacent triangles
		if(links[i].t1 == -1 || links[i].t2 == -1)
			continue;

		t1 = triangles_[links[i].t1];
		t2 = triangles_[links[i].t2];
		i0 = links[i].i0;
		i1 = links[i].i1;

		neighbors[links[i].t1].push_back(links[i].t2);
		neighbors[links[i].t2].push_back(links[i].t1);

		// Permute the indices of the first triangle so that the 0th and 2nd index are the shared edge
		while( ! ((t1.i0 == i0 && t1.i2 == i1) || (t1.i0 == i1 && t1.i2 == i0)))
			detail::permute_triangle_vertices(t1);

		// Permute the indices of the second triangle so that the 0th and 1st index are the shared edge
		while( ! ((t2.i0 == i0 && t2.i1 == i1) || (t2.i0 == i1 && t2.i1 == i0)))
			detail::permute_triangle_vertices(t2);

		tri_pair.i0 = t1.i0;
		tri_pair.i1 = t1.i1;
		tri_pair.i2 = t1.i2;
		tri_pair.i3 = t2.i2;

		// Compute sine and cosine of the angle between the normals
		tri::cos_sin_of_angle_between_pair(
				vertices_[tri_pair.i0],
				vertices_[tri_pair.i1],
				vertices_[tri_pair.i2],
				vertices_[tri_pair.i3],
				tri_pair.cos_theta0,
				tri_pair.sin_theta0);

		triangle_pairs_.push_back(tri_pair);
	}

	// Create triangle quadruples (triangle+neighbors)
	for(typename std::map<plint, std::vector<plint> >::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
		TriangleQuadruple quad;

		// Indices 0, 1, 2 correspond to the main triangle
		Triangle<T> & main_tri = triangles_[it->first];
		quad.i0 = main_tri.i0;
		quad.i1 = main_tri.i1;
		quad.i2 = main_tri.i2;

		for(typename std::vector<plint>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
			Triangle<T> & t = triangles_[*it2];

			// Check which edge the neighboring triangle shares with the main triangle
			plint edge;
			if((edge = detail::find_missing_vertex(main_tri.i0, main_tri.i1, t)) != -1)
				quad.i5 = edge;
			if((edge = detail::find_missing_vertex(main_tri.i1, main_tri.i2, t)) != -1)
				quad.i4 = edge;
			if((edge = detail::find_missing_vertex(main_tri.i2, main_tri.i0, t)) != -1)
				quad.i3 = edge;
		}
	}
}

template<class T>
ParticleShape<T>::~ParticleShape()
{

}

template<class T>
void ParticleShape<T>::print_statistics() const
{
	std::cout << "Mesh statistics:" << std::endl;
	std::cout << "  #vertices: " << vertices_.size() << std::endl;

	T A_avg = T();
	T A_std = T();
	for(triangle_const_iterator it = triangles_.begin(); it != triangles_.end(); ++it)
		A_avg += it->area;
	A_avg /= triangles_.size();

	for(triangle_const_iterator it = triangles_.begin(); it != triangles_.end(); ++it)
		A_std += util::sqr(it->area - A_avg);
	A_std = std::sqrt(A_std / triangles_.size());

	std::cout << "  #elements: " << triangles_.size() << ", average area: " << A_avg << ", std: " << A_std << std::endl;

	T l_avg = T();
	T l_std = T();
	for(link_const_iterator it = links.begin(); it != links.end(); ++it)
		l_avg += it->length;
	l_avg /= links.size();

	for(link_const_iterator it = links.begin(); it != links.end(); ++it)
		l_std += util::sqr(it->length - l_avg);
	l_std = std::sqrt(l_std / links.size());

	std::cout << "  #links: " << links.size() << ", average length: " << l_avg << ", std: " << l_std << std::endl;

	std::cout << "  Surface area: " << area << ", volume: " << volume << std::endl;
}

template <class T>
void ParticleShape<T>::compute_properties()
{
    const T oneDiv6 = (T)(1.0/6.0);
    const T oneDiv24 = (T)(1.0/24.0);
    const T oneDiv60 = (T)(1.0/60.0);
    const T oneDiv120 = (T)(1.0/120.0);

    // order:  1, x, y, z, x^2, y^2, z^2, xy, yz, zx
    T integral[10] = { (T)0.0, (T)0.0, (T)0.0, (T)0.0,
        (T)0.0, (T)0.0, (T)0.0, (T)0.0, (T)0.0, (T)0.0 };

    for(triangle_const_iterator it = triangles_begin(); it != triangles_end(); ++it) {
        // Get vertices of triangle i.
        const Array<T, 3> & v0 = get_vertex(it->i0);
        const Array<T, 3> & v1 = get_vertex(it->i1);
        const Array<T, 3> & v2 = get_vertex(it->i2);

        // Get cross product of edges and normal vector.
        Array<T, 3> V1mV0 = v1 - v0;
        Array<T, 3> V2mV0 = v2 - v0;
        Array<T, 3> N = crossProduct(V1mV0, V2mV0);

        // Store normal
        //T normal_norm = norm(N);
        //nodes[i].normal = N / normal_norm;

        // Compute area
        //nodes[i].area = 0.5 * normal_norm;

        // Compute integral terms.
        T tmp0, tmp1, tmp2;
        T f1x, f2x, f3x, g0x, g1x, g2x;
        tmp0 = v0[0] + v1[0];
        f1x = tmp0 + v2[0];
        tmp1 = v0[0]*v0[0];
        tmp2 = tmp1 + v1[0]*tmp0;
        f2x = tmp2 + v2[0]*f1x;
        f3x = v0[0]*tmp1 + v1[0]*tmp2 + v2[0]*f2x;
        g0x = f2x + v0[0]*(f1x + v0[0]);
        g1x = f2x + v1[0]*(f1x + v1[0]);
        g2x = f2x + v2[0]*(f1x + v2[0]);

        T f1y, f2y, f3y, g0y, g1y, g2y;
        tmp0 = v0[1] + v1[1];
        f1y = tmp0 + v2[1];
        tmp1 = v0[1]*v0[1];
        tmp2 = tmp1 + v1[1]*tmp0;
        f2y = tmp2 + v2[1]*f1y;
        f3y = v0[1]*tmp1 + v1[1]*tmp2 + v2[1]*f2y;
        g0y = f2y + v0[1]*(f1y + v0[1]);
        g1y = f2y + v1[1]*(f1y + v1[1]);
        g2y = f2y + v2[1]*(f1y + v2[1]);

        T f1z, f2z, f3z, g0z, g1z, g2z;
        tmp0 = v0[2] + v1[2];
        f1z = tmp0 + v2[2];
        tmp1 = v0[2]*v0[2];
        tmp2 = tmp1 + v1[2]*tmp0;
        f2z = tmp2 + v2[2]*f1z;
        f3z = v0[2]*tmp1 + v1[2]*tmp2 + v2[2]*f2z;
        g0z = f2z + v0[2]*(f1z + v0[2]);
        g1z = f2z + v1[2]*(f1z + v1[2]);
        g2z = f2z + v2[2]*(f1z + v2[2]);

        // Update integrals.
        integral[0] += N[0]*f1x;
        integral[1] += N[0]*f2x;
        integral[2] += N[1]*f2y;
        integral[3] += N[2]*f2z;
        integral[4] += N[0]*f3x;
        integral[5] += N[1]*f3y;
        integral[6] += N[2]*f3z;
        integral[7] += N[0]*(v0[1]*g0x + v1[1]*g1x + v2[1]*g2x);
        integral[8] += N[1]*(v0[2]*g0y + v1[2]*g1y + v2[2]*g2y);
        integral[9] += N[2]*(v0[0]*g0z + v1[0]*g1z + v2[0]*g2z);
    }

    integral[0] *= oneDiv6;
    integral[1] *= oneDiv24;
    integral[2] *= oneDiv24;
    integral[3] *= oneDiv24;
    integral[4] *= oneDiv60;
    integral[5] *= oneDiv60;
    integral[6] *= oneDiv60;
    integral[7] *= oneDiv120;
    integral[8] *= oneDiv120;
    integral[9] *= oneDiv120;

    // mass
    //volume = integral[0];

    // center of mass
    //center = Array<T, 3>(integral[1], integral[2], integral[3])/volume;

    // inertia relative to center of mass
    moment_of_inertia(0, 0) = integral[5] + integral[6];
    moment_of_inertia(0, 1) = -integral[7];
    moment_of_inertia(0, 2) = -integral[9];
    moment_of_inertia(1, 0) = moment_of_inertia(0, 1);
    moment_of_inertia(1, 1) = integral[4] + integral[6];
    moment_of_inertia(1, 2) = -integral[8];
    moment_of_inertia(2, 0) = moment_of_inertia(0, 2);
    moment_of_inertia(2, 1) = moment_of_inertia(1, 2);
    moment_of_inertia(2, 2) = integral[4] + integral[5];
}

template <class T>
void ParticleShape<T>::compute_bounding_volumes()
{
    // Minimal bounding sphere radius
    radius = 0;

    // Minimal axis aligned bounding box
    bounding_box.x0 = bounding_box.y0 = bounding_box.z0 = std::numeric_limits<T>::max();
    bounding_box.x1 = bounding_box.y1 = bounding_box.z1 = -std::numeric_limits<T>::max();

	for(unsigned int i = 0; i < vertices_.size(); ++i) {
		T dist_sqr = util::sqr(vertices_[i][0]) + util::sqr(vertices_[i][1]) + util::sqr(vertices_[i][2]);
		if(dist_sqr > radius)
			radius = dist_sqr;

		bounding_box.x0 = std::min(bounding_box.x0, vertices_[i][0]);
		bounding_box.y0 = std::min(bounding_box.y0, vertices_[i][1]);
		bounding_box.z0 = std::min(bounding_box.z0, vertices_[i][2]);
		bounding_box.x1 = std::max(bounding_box.x1, vertices_[i][0]);
		bounding_box.y1 = std::max(bounding_box.y1, vertices_[i][1]);
		bounding_box.z1 = std::max(bounding_box.z1, vertices_[i][2]);
	}

	radius = std::sqrt(radius);

	// Save the integer counterpart for the bounding box
	bounding_box_i.x0 = std::floor(bounding_box.x0);
	bounding_box_i.x1 = std::floor(bounding_box.x1) + 1;
	bounding_box_i.y0 = std::floor(bounding_box.y0);
	bounding_box_i.y1 = std::floor(bounding_box.y1) + 1;
	bounding_box_i.z0 = std::floor(bounding_box.z0);
	bounding_box_i.z1 = std::floor(bounding_box.z1) + 1;
}



} /* namespace fsi */

} /* namespace plb */

#endif /* PARTICLESHAPE_HH_ */
