/*
 * Periodicity.h
 *
 *  Created on: Apr 8, 2014
 *      Author: niber
 */

#ifndef PERIODICITY_H_
#define PERIODICITY_H_

namespace plb {

namespace fsi {

/*
 * These classes make it possible to evaluate vector operations in modular arithmetic.
 * The Periodicity and Arithmetic3D classes are templatized with 3 booleans corresponding
 * to the periodicity of the domain. At compile time these are known and we can thus generate
 * the minimal amount of code for periodicity imposing.
 */

namespace detail {

// Minimal modulus (e.g. -3%7 === -3, and not 4). The result always gives the
// smallest possible absolute value
inline float min_mod(float x, float N, float N_inv)
{
	return x - floorf(x*N_inv + 0.5)*N;
}

inline double min_mod(double x, double N, double N_inv)
{
	return x - floor(x*N_inv + 0.5)*N;
}

// Positive modulus (e.g. -3%7 === 4, and not -3). The result is always positive.
inline double pos_mod(float x, float N, float N_inv)
{
	return x = x - floorf(x*N_inv)*N;
}

inline double pos_mod(double x, double N, double N_inv)
{
	return x = x - floor(x*N_inv)*N;
}

// Forward declaration
template<class T, bool val>
struct Arithmetic1D
{
	static void remap_pos(T & x, const T & x0, const T & N, const T & N_inv);
	static T min_mod_diff(const T & x, const T & y, const T & N, const T & N_inv);
	static void remap_int(plint & i, const plint & offset, plint N);
	static void shift_toward_value(const T & ref, T & x, const T & N, const T & N_inv);
};

// Specialization for a periodic direction
template<class T>
struct Arithmetic1D<T, true> {
	// Remap position from [-infty, infty] to [x0, x0+N)
	static void remap_pos(T & x, const T & x0, const T & N, const T & N_inv)
	{
		x = x - std::floor((x-x0)*N_inv)*N;
	}

	// Minimal difference in modular sense
	static T min_mod_diff(const T & x, const T & y, const T & N, const T & N_inv)
	{
		return min_mod(x-y, N, N_inv);
	}

	// Remap an int from (-inf, inf) to [0, N-1]
	static void remap_int(plint & i, const plint & offset, plint N)
	{
		i = (N + ((i+offset)%N)) % N - offset;
	}

	// Add/subtract a multiple of N to x so that it is as close as possible to ref (the result might be both positive and negative)
	static void shift_toward_value(const T & ref, T & x, const T & N, const T & N_inv)
	{
		x = x - std::floor((x-ref)*N_inv + 0.5)*N;
	}
};

// Non-periodic specialization
// The methods do not perform any operations and the compiler should optimize away calls.
template<class T>
struct Arithmetic1D<T, false> {
	static void remap_pos(T & x, const T & x0, const T & N, const T & N_inv) {  }
	static T min_mod_diff(const T & x, const T & y, const T & N, const T & N_inv) { return x - y; }
	static void remap_int(plint & i, const plint & offset, plint N) { }
	static void shift_toward_value(const T & ref, T & x, const T & N, const T & N_inv) { }
};

template<class T> T pow2(const T & x) { return x*x; }

}

template<class T, bool pX, bool pY, bool pZ>
class Arithmetic3D {
public:
	Arithmetic3D(const Box3D & dom) :
	domain(dom),
	nxi(dom.getNx()), nyi(dom.getNy()), nzi(dom.getNz()),
	x0(dom.x0-0.5), y0(dom.y0-0.5), z0(dom.z0-0.5)
	{
		nx = (T) (nxi);
		ny = (T) (nyi);
		nz = (T) (nzi);
		nx_inv = (T) 1 / nx;
		ny_inv = (T) 1 / ny;
		nz_inv = (T) 1 / nz;
	}

	// Remap an int to [0, nx-1], [0, ny-1], [0, nz-1], respectively.
	// If a direction is non-periodic, i is unchanged.
	void remap_index_x(plint & i, const plint & offset) const { detail::Arithmetic1D<T, pX>::remap_int(i, offset, nxi); }
	void remap_index_y(plint & i, const plint & offset) const { detail::Arithmetic1D<T, pY>::remap_int(i, offset, nyi); }
	void remap_index_z(plint & i, const plint & offset) const { detail::Arithmetic1D<T, pZ>::remap_int(i, offset, nzi); }

	// Remap a position to [x0, x0+nx) x [y0, y0+ny) x [z0, z0+nz)
	// A component of the array remains unchanged if the corresponding direction is non-periodic.
	void remap_position(Array<T, 3> & pos) const
	{
		detail::Arithmetic1D<T, pX>::remap_pos(pos[0], x0, nx, nx_inv);
		detail::Arithmetic1D<T, pY>::remap_pos(pos[1], y0, ny, ny_inv);
		detail::Arithmetic1D<T, pZ>::remap_pos(pos[2], z0, nz, nz_inv);
	}

	// Push the array pos as close to ref as possible, only adding/subtracting multiples of nx, ny, nz
	// to each component of x, respectively.
	void shift_periodically_to_minimize_distance_to(const Array<T, 3> & ref, Array<T, 3> & pos) const
	{
		detail::Arithmetic1D<T, pX>::shift_toward_value(ref[0], pos[0], nx, nx_inv);
		detail::Arithmetic1D<T, pY>::shift_toward_value(ref[1], pos[1], ny, ny_inv);
		detail::Arithmetic1D<T, pZ>::shift_toward_value(ref[2], pos[2], nz, nz_inv);
	}

	// Minimal distance between two values (in modular sense) in a given direction
	T dist_x(const T & p0, const T & p1) const { return detail::Arithmetic1D<T, pX>::min_mod_diff(p0, p1, nx, nx_inv); }
	T dist_y(const T & p0, const T & p1) const { return detail::Arithmetic1D<T, pY>::min_mod_diff(p0, p1, ny, ny_inv); }
	T dist_z(const T & p0, const T & p1) const { return detail::Arithmetic1D<T, pZ>::min_mod_diff(p0, p1, nz, nz_inv); }

	// Distance between two 2D vectors
	T dist_sqr_xy(const Array<T, 2> & v1, const Array<T, 2> & v2) const
	{
		return detail::pow2(detail::Arithmetic1D<T, pX>::min_mod_diff(v1[0], v2[0], nx, nx_inv)) +
				detail::pow2(detail::Arithmetic1D<T, pY>::min_mod_diff(v1[1], v2[1], ny, ny_inv));
	}

	T dist_sqr_xz(const Array<T, 2> & v1, const Array<T, 2> & v2) const
	{
		return detail::pow2(detail::Arithmetic1D<T, pX>::min_mod_diff(v1[0], v2[0], nx, nx_inv)) +
				detail::pow2(detail::Arithmetic1D<T, pZ>::min_mod_diff(v1[1], v2[1], nz, nz_inv));
	}

	T dist_sqr_yz(const Array<T, 2> & v1, const Array<T, 2> & v2) const
	{
		return detail::pow2(detail::Arithmetic1D<T, pY>::min_mod_diff(v1[0], v2[0], ny, ny_inv)) +
				detail::pow2(detail::Arithmetic1D<T, pZ>::min_mod_diff(v1[1], v2[1], nz, nz_inv));
	}

	// Distance between two 3D vectors
	T dist_sqr(const Array<T, 3> & v1, const Array<T, 3> & v2) const
	{
		return detail::pow2(detail::Arithmetic1D<T, pX>::min_mod_diff(v1[0], v2[0], nx, nx_inv)) +
				detail::pow2(detail::Arithmetic1D<T, pY>::min_mod_diff(v1[1], v2[1], ny, ny_inv)) +
				detail::pow2(detail::Arithmetic1D<T, pZ>::min_mod_diff(v1[2], v2[2], nz, nz_inv));
	}

	T dist(const Array<T, 3> & v1, const Array<T, 3> & v2) const
	{
		return std::sqrt(dist_sqr(v1, v2));
	}

	Array<T, 3> vec_diff(const Array<T, 3> & v1, const Array<T, 3> & v2) const
	{
		return Array<T, 3>(
				detail::Arithmetic1D<T, pX>::min_mod_diff(v1[0], v2[0], nx, nx_inv),
				detail::Arithmetic1D<T, pY>::min_mod_diff(v1[1], v2[1], ny, ny_inv),
				detail::Arithmetic1D<T, pZ>::min_mod_diff(v1[2], v2[2], nz, nz_inv)
		);
	}

private:
	pluint nxi, nyi, nzi;
	T x0, y0, z0;
	T nx, ny, nz, nx_inv, ny_inv, nz_inv;
	Box3D domain;
};

template<class T>
struct NormalArithmetic : public Arithmetic3D<T, false, false, false> {
	NormalArithmetic() : Arithmetic3D<T, false, false, false>(Box3D(0, 0, 0, 0, 0, 0)) { }
};

template<class T, bool pX, bool pY, bool pZ>
class Periodicity3D {
public:
	typedef Arithmetic3D<T, pX, pY, pZ> ArithmeticType;

	static bool get_x() { return pX; }
	static bool get_y() { return pY; }
	static bool get_z() { return pZ; }

	static ArithmeticType create_arithmetic(const Box3D & domain)
	{
		return Arithmetic3D<T, pX, pY, pZ>(domain);
	}

private:

};

}

}




#endif /* PERIODICITY_H_ */
