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

namespace plb {

namespace fsi {

template<class T>
ParticleShape<T>::ParticleShape(Array<T, 3> * vertices, unsigned int * indices, unsigned int num_vertices, unsigned int num_elems)
: vertices_(vertices), indices(indices), num_triangles(num_elems), num_vertices(num_vertices), nodes(num_elems), scale(1.0)
{
	// Compute mass properties
	compute_properties();

	// Find the neighbors of each triangle
	find_mesh_neighbors();

	// Compute bounding volumes (box and sphere)
	compute_bounding_volumes();

	hash_triangles();
}

template<class T>
void ParticleShape<T>::find_mesh_neighbors()
{
	const unsigned int nobf = 3;

	// This code is written by Minh Duon-Quang as a part of the Slilab library
	int *ia = (int*) malloc( (num_triangles+1)*sizeof(int) );
	int *ja = (int*) malloc( (num_triangles*nobf)*sizeof(int) );

	int *nptr, *nind, *mark, *jat, *iat;
	int i, m, el, ii, j, jj, k, kk, kkk, node, nedges, start, end, nz, nzmax, nzmaxa;
	int mask = (1<<11)-1;
	int ind[200], wgt[200];
	int mgcnum = 2; // 2D surface only

	mark = (int*)malloc((mask+1)*sizeof(int));
	for(i=0; i<=mask; i++)
		mark[i] = -1;

	/* Construct the node-element list first */
	nptr = (int*)malloc((num_vertices+1)*sizeof(int));

	for(i=0;i<=num_vertices;i++) nptr[i]=0;

	for (el=0; el<num_triangles; el++)
		for(i=0; i<3; i++)
			nptr[indices[el*nobf+i]]++;

	for (i=1; i < num_vertices; i++)
		nptr[i] += nptr[i-1];
	for (i=num_vertices; i > 0; i--)
		nptr[i] = nptr[i-1];
	nptr[0] = 0;

	nind = (int*)malloc(nptr[num_vertices]*sizeof(int));
	for (i=0; i<num_triangles; i++) {
		for (j=0; j<nobf; j++)
			nind[nptr[indices[i*nobf+j]]++] = i;
	}
	for (i=num_vertices; i>0; i--)
		nptr[i] = nptr[i-1];
	nptr[0] = 0;

	for (el=0; el<num_triangles; el++)
		ia[el] = nobf*el;

	for (el=0; el<num_triangles; el++) {
		for(m=i=0; i<nobf; i++) {
			node = indices[el*nobf+i];
			for(k=nptr[node+1]-1; k>=nptr[node]; k--) {
				if((kk = nind[k]) <= el)
					break;

				kkk = kk & mask;

				if((j=mark[kkk]) == -1) {
					ind[m] = kk;
					wgt[m] = 1;
					mark[kkk] = m++;
				}
				else if (ind[j] == kk) {
					wgt[j]++;
				}
				else {
					for(jj=0; jj<m; jj++) {
						if(ind[jj] == kk){
							wgt[jj]++;
							break;
						}
					}
					if(jj == m) {
						ind[m] = kk;
						wgt[m++] = 1;
					}
				}
			}
		}

		for(i=0; i<m; i++) {
			if(wgt[i] == mgcnum) {
				k = ind[i];
				ja[ia[el]++] = k;
				ja[ia[ k]++] = el;
			}
			mark[ind[i]&mask] = -1;
		}
	}

	for(j=i=ii=0; i < num_triangles; i++) {
		for(k=nobf*i; k<ia[i]; k++, j++)
			ja[j] = ja[k];
		ia[ii++] = j;
	}
	for(i=ii; i>0; i--)
		ia[i] = ia[i-1];
	ia[0] = 0;

	free(mark);
	free(nptr);
	free(nind);

	for(int i=0; i<num_triangles; i++) {
		T lambda = T();
		this->nodes[i].neighbors.clear();
		for (int idx=0, j=ia[i]; j<ia[i+1]; j++) {
			Array<T,3> vds = this->nodes[i].centroid - this->nodes[ja[j]].centroid;
			Array<T,3> vdT = this->nodes[i].normal - this->nodes[ja[j]].normal;
			lambda += norm(vdT) / norm(vds);
			this->nodes[i].neighbors.push_back(ja[j]);
		}
		lambda /= (ia[i+1]-ia[i]);
		this->nodes[i].lambda = lambda;
	}
	free(ia);
	free(ja);
}

template<class T>
ParticleShape<T>::~ParticleShape()
{
	delete [] vertices_;
	delete [] indices;
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

    const unsigned int* index = indices;
    for(plint i = 0; i < num_triangles; i++)
    {
    	nodes[i].id = i;

        // Get vertices of triangle i.
    	nodes[i].vertices_[0] = vertices_[*index++];
    	nodes[i].vertices_[1] = vertices_[*index++];
    	nodes[i].vertices_[2] = vertices_[*index++];

        Array<T, 3> v0 = scale * nodes[i].vertices_[0];
        Array<T, 3> v1 = scale * nodes[i].vertices_[1];
        Array<T, 3> v2 = scale * nodes[i].vertices_[2];

        // Compute centroid
        nodes[i].centroid = (v0 + v1 + v2) / (T) 3.0;

        // Get cross product of edges and normal vector.
        Array<T, 3> V1mV0 = v1 - v0;
        Array<T, 3> V2mV0 = v2 - v0;
        Array<T, 3> N = crossProduct(V1mV0, V2mV0);

        // Store normal
        T normal_norm = norm(N);
        nodes[i].normal = N / normal_norm;

        // Compute area
        nodes[i].area = 0.5 * normal_norm;

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
    volume = integral[0];

    // center of mass
    center = Array<T, 3>(integral[1], integral[2], integral[3])/volume;

    // inertia relative to center of mass
    moment_of_inertia(0, 0) = integral[5] + integral[6] /*-
    		volume*(center[1]*center[1] + center[2]*center[2])*/;
    moment_of_inertia(0, 1) = -integral[7] /*+ volume*center[0]*center[1]*/;
    moment_of_inertia(0, 2) = -integral[9] /*+ volume*center[2]*center[0]*/;
    moment_of_inertia(1, 0) = moment_of_inertia(0, 1);
    moment_of_inertia(1, 1) = integral[4] + integral[6] /*-
    		volume*(center[2]*center[2] + center[0]*center[0])*/;
    moment_of_inertia(1, 2) = -integral[8] /*+ volume*center[1]*center[2]*/;
    moment_of_inertia(2, 0) = moment_of_inertia(0, 2);
    moment_of_inertia(2, 1) = moment_of_inertia(1, 2);
    moment_of_inertia(2, 2) = integral[4] + integral[5] /*-
    		volume*(center[0]*center[0] + center[1]*center[1])*/;
}

template <class T>
void ParticleShape<T>::compute_bounding_volumes()
{
    // Minimal bounding sphere radius
    radius = 0;

    // Minimal axis aligned bounding box
    bounding_box.x0 = bounding_box.y0 = bounding_box.z0 = std::numeric_limits<T>::max();
    bounding_box.x1 = bounding_box.y1 = bounding_box.z1 = -std::numeric_limits<T>::max();

	for(unsigned int i = 0; i < num_vertices; ++i) {
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

template <class T>
void ParticleShape<T>::hash_triangles()
{
	// Puts all triangles into a grid for efficient ray tracing and nearest neighbor search
	// It should however be noted that the current algorithm is quite conservative and could be optimized

	hash_dim.x = bounding_box_i.x1 - bounding_box_i.x0;
	hash_dim.y = bounding_box_i.y1 - bounding_box_i.y0;
	hash_dim.z = bounding_box_i.z1 - bounding_box_i.z0;

	// Allocate grid
	triangle_hash.resize(hash_dim.x * hash_dim.y * hash_dim.z);

	for(plint i = 0; i < num_triangles; ++i) {
		triangle_hash[hash_position(nodes[i].centroid)].push_back(&(nodes[i]));
	}
}

template <class T>
plint ParticleShape<T>::hash_position(const Array<T, 3> & pos) const
{
	return (std::floor(pos[0]-bounding_box.x0)*hash_dim.y +
			 std::floor(pos[1]-bounding_box.y0))*hash_dim.z +
			 std::floor(pos[2]-bounding_box.z0);
}

template <class T>
plint ParticleShape<T>::get_closest_element_facing_in_direction(const Array<T, 3> & r, const Array<T, 3> & dir, T cutoff) const
{
	T dist_sqr;
	T min_dist = std::numeric_limits<T>::max();
	plint ret = -1;

	if(cutoff > 0.5*radius) {
		// Large cutoff radius, fall back to O(N) test
		for(unsigned int i = 0; i < nodes.size(); ++i) {
			if(dot(dir, nodes[i].normal) > 0) {
				dist_sqr = normSqr(r - nodes[i].centroid);
				if(dist_sqr < min_dist) {
					ret = i;
					min_dist = dist_sqr;
				}
			}
		}
	} else {
		// Find bounding box
		plint x_min = std::max((plint)(r[0] - cutoff - bounding_box.x0), (plint)0);
		plint x_max = std::min((plint)(r[0] + cutoff - bounding_box.x0), hash_dim.x-1);
		plint y_min = std::max((plint)(r[1] - cutoff - bounding_box.y0), (plint)0);
		plint y_max = std::min((plint)(r[1] + cutoff - bounding_box.y0), hash_dim.y-1);
		plint z_min = std::max((plint)(r[2] - cutoff - bounding_box.z0), (plint)0);
		plint z_max = std::min((plint)(r[2] + cutoff - bounding_box.z0), hash_dim.z-1);

		for(plint i = x_min; i <= x_max; ++i)
			for(plint j = y_min; j <= y_max; ++j)
				for(plint k = z_min; k <= z_max; ++k) {
					const std::vector<ParticleShapeNode<T> *> & ns = triangle_hash[(i*hash_dim.y + j)*hash_dim.z + k];
					for(plint l = 0; l < ns.size(); ++l) {
						if(dot(dir, ns[l]->normal) > 0) {
							dist_sqr = normSqr(r - ns[l]->centroid);
							if(dist_sqr < min_dist) {
								ret = ns[l]->id;
								min_dist = dist_sqr;
							}
						}
					}
				}
	}

	return ret;
}

template <class T>
void ParticleShape<T>::trace_ray(const Array<T, 3> & pos, const Array<T, 3> & dir, plint & id, T & dist) const
{
	T epsilon = 0.00001;
	dist = std::numeric_limits<T>::max();
	id = -1;

	for(unsigned int i = 0; i < nodes.size(); ++i) {
		if(dot(dir, nodes[i].normal) > 0)
			continue;

		const Array<T, 3> e1 = nodes[i].vertices_[1]-nodes[i].vertices_[0];
		const Array<T, 3> e2 = nodes[i].vertices_[2]-nodes[i].vertices_[0];
		const Array<T, 3> q  = crossProduct(dir,e2);
		const T a  = dot(e1,q);

		// Check if the vector is parallel to the plane (the intersection is at infinity)
		if (a>-epsilon && a<epsilon)
			continue;

		const T f = 1/a;
		const Array<T, 3> s = pos - nodes[i].vertices_[0];
		const T u = f*dot(s,q);

		// Check if the intersection is outside of the triangle
		if (u<0.0)
			continue;

		const Array<T, 3> r = crossProduct(s,e1);
		const T v = f*dot(dir,r);

		// Check if the intersection is outside of the triangle
		if(v<0.0 || u+v>1.0)
			continue;

		const T t = f*dot(e2,r);
		if(t < dist) {
			dist = t;
			id = i;
		}
	}
}

template <class T>
Array<T, 3> ParticleShape<T>::get_average_distance_within_sphere(const Array<T, 3> & v, T r) const
{
	T dist_sqr;
	T r_sqr = r*r;
	T A_tot = 0;
	Array<T, 3> r_tot; r_tot.resetToZero();

	for(unsigned int i = 0; i < nodes.size(); ++i) {
		Array<T, 3> dr = nodes[i].centroid - v;
		if(normSqr(dr) < r_sqr) {
			r_tot += dr * nodes[i].area;
			A_tot += nodes[i].area;
		}
	}

	if(A_tot < 1e-3)
		return Array<T, 3>(1000000, 1000000, 1000000);

	r_tot /= A_tot;
	return r_tot;
}

} /* namespace fsi */

} /* namespace plb */

#endif /* PARTICLESHAPE_HH_ */
