/*
 * linalg.h
 *
 *  Created on: Jun 25, 2015
 *      Author: niber
 */

#ifndef LINALG_H_
#define LINALG_H_
#include "core/array.h"
#include "Quaternion.h"
#include "Matrix3.h"
#include <cmath>

namespace plb {

namespace fsi {

namespace linalg {

template<class T>
void diagonalize(Matrix<T, 3> & A, Matrix<T, 3> & Q, Array<T, 3> & w)
{
	// ----------------------------------------------------------------------------
	// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
	// matrix A using the Jacobi algorithm.
	// The upper triangular part of A is destroyed during the calculation,
	// the diagonal elements are read but not destroyed, and the lower
	// triangular elements are not referenced at all.
	// ----------------------------------------------------------------------------
	// Parameters:
	//   A: The symmetric input matrix
	//   Q: Storage buffer for eigenvectors
	//   w: Storage buffer for eigenvalues
	// ----------------------------------------------------------------------------
	// Return value:
	//   0: Success
	//  -1: Error (no convergence)
	// ----------------------------------------------------------------------------
	const int n = 3;
	T sd, so;                  // Sums of diagonal resp. off-diagonal elements
	T s, c, t;                 // sin(phi), cos(phi), tan(phi) and temporary storage
	T g, h, z, theta;          // More temporary storage
	T thresh;

	// Initialize Q to the identitity matrix
	for (int i=0; i < n; i++) {
		Q(i, i) = 1.0;
		for (int j=0; j < i; j++)
			Q(i ,j) = Q(j, i) = 0.0;
	}

	// Initialize w to diag(A)
	for (int i=0; i < n; i++)
		w[i] = A(i, i);

	// Calculate SQR(tr(A))
	sd = 0.0;
	for (int i=0; i < n; i++)
		sd += std::abs(w[i]);
	sd = sd * sd;

	// Main iteration loop
	for (int nIter=0; nIter < 50; nIter++) {
		// Test for convergence
		so = 0.0;
		for (int p=0; p < n; p++)
			for (int q=p+1; q < n; q++)
				so += std::abs(A(p, q));

			if (so == 0.0)
				return;

			if (nIter < 4)
				thresh = 0.2 * so / util::sqr(n);
			else
				thresh = 0.0;

			// Do sweep
			for (int p=0; p < n; p++)
				for (int q=p+1; q < n; q++) {
					g = 100.0 * std::abs(A(p, q));
					if (nIter > 4  &&  std::abs(w[p]) + g == std::abs(w[p])
						   &&  std::abs(w[q]) + g == std::abs(w[q])) {
						A(p, q) = 0.0;
					} else if (std::abs(A(p, q)) > thresh) {
						// Calculate Jacobi transformation
						h = w[q] - w[p];
						if (std::abs(h) + g == std::abs(h)) {
							t = A(p, q) / h;
						} else {
							theta = 0.5 * h / A(p, q);
							if (theta < 0.0)
								t = -1.0 / (std::sqrt(1.0 + util::sqr(theta)) - theta);
							else
								t = 1.0 / (std::sqrt(1.0 + util::sqr(theta)) + theta);
						}
						c = 1.0/std::sqrt(1.0 + util::sqr(t));
						s = t * c;
						z = t * A(p, q);

						// Apply Jacobi transformation
						A(p, q) = 0.0;
						w[p] -= z;
						w[q] += z;
						for (int r=0; r < p; r++) {
							t = A(r, p);
							A(r, p) = c*t - s*A(r, q);
							A(r, q) = s*t + c*A(r, q);
						}

						for (int r=p+1; r < q; r++) {
							t = A(p, r);
							A(p, r) = c*t - s*A(r, q);
							A(r, q) = s*t + c*A(r, q);
						}

						for (int r=q+1; r < n; r++) {
							t = A(p, r);
							A(p, r) = c*t - s*A(q, r);
							A(q, r) = s*t + c*A(q, r);
						}

						// Update eigenvectors
						for (int r=0; r < n; r++) {
							t = Q(r, p);
							Q(r, p) = c*t - s*Q(r, q);
							Q(r, q) = s*t + c*Q(r, q);
						}
					}
				}
	}
}


} /* namespace linalg */


} /* namespace fsi */

} /* namespace plb */



#endif /* LINALG_H_ */
