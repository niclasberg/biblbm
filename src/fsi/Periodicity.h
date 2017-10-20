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
// May suffer from round off errors
inline float min_mod(float x, float N, float N_inv)
{
	return x - floorf(x*N_inv + 0.5)*N;
}

inline double min_mod(double x, double N, double N_inv)
{
	return x - std::floor(x*N_inv + 0.5)*N;
}

// Positive modulus (e.g. -3%7 === 4, and not -3). The result is always positive.
//
template<class T>
inline T pos_mod(T x, T N)
{
	return (T) std::fmod(x + N, N);
}

// Default arithmetic, does nothing
template<class T, bool val>
class Arithmetic1D
{
public:
	Arithmetic1D(plint x0, plint x1)
	: x0_(x0), x1_(x1), N_(x1 - x0 + 1)
	{ }

	plint get_N() const { return N_; }
	void remap_pos(T & x) const { }
	plint int_diff(const plint & i0, const plint & i1) const { return i0-i1; }
	T min_mod_diff(const T & x, const T & y) const { return x-y; }
	plint remap_int(const plint & i) const { return i; }
	plint remap_int(const plint & i, const plint & offset) const { return i; }
	void shift_toward_value(const T & ref, T & x) const { }

private:
	plint x0_, x1_, N_;
};

// Specialization for a periodic direction
template<class T>
class Arithmetic1D<T, true> {
public:
	Arithmetic1D(plint x0, plint x1)
	: i0_(x0), i1_(x1), N_(x1 - x0 + 1)
	{
		// Domain in floating point representation
		x0_ = (T)i0_ - 0.5;
		x1_ = (T)i1_ + 0.5;
		L_ = (T) N_;
		Linv_ = 1. / L_;
	}

	plint get_N() const { return N_; }

	// Remap position from [-infty, infty] to [i0-0.5, i1+0.5)
	void remap_pos(T & x) const
	{
		x = x - L_*std::floor((x-x0_)*Linv_);
		//x = x0_ + std::fmod(x-x0_ + 2.*L_, L_);
	}

	// Integer difference (in interval [0, N_-1])
	plint int_diff(const plint & i0, const plint & i1) const
	{
		return ((i0-i1 + 100*N_)%N_);
	}

	// Minimal difference in modular sense
	T min_mod_diff(const T & x, const T & y) const
	{
		return min_mod(x-y, L_, Linv_);
	}

	// Remap an int from (-inf, inf) to [i0, N+i0-1]
	plint remap_int(plint i) const
	{
		if (i < i0_)
			i += N_ * ((i0_ - i) / N_ + 1);
		return i0_ + (i - i0_) % N_;
	}

	// Remap (i-offset) to [i0, N+i0-1]
	plint remap_int(const plint & i, const plint & offset) const
	{
		plint val = ((i+offset-i0_ + 100*N_) % N_) - offset + i0_;
		return (val < i0_) ? val+N_ : val;
	}

	// Add/subtract a multiple of N to x so that it is as close as possible to ref (the result might be both positive and negative)
	void shift_toward_value(const T & ref, T & x) const
	{
		x = x - std::floor((x-ref)*Linv_ + 0.5)*L_;
	}
private:
	plint i0_, i1_, N_;
	T x0_, x1_, L_, Linv_;
};

template<class T> T pow2(const T & x) { return x*x; }

}

template<class T, bool pX, bool pY, bool pZ>
class Arithmetic3D {
public:
	Arithmetic3D(const Box3D & dom) :
	domain(dom),
	arithmetic_x(dom.x0, dom.x1),
	arithmetic_y(dom.y0, dom.y1),
	arithmetic_z(dom.z0, dom.z1)
	{

	}

	// Getters
	bool periodic_x() const { return pX; }
	bool periodic_y() const { return pY; }
	bool periodic_z() const { return pZ; }
	plint get_nx() const { return arithmetic_x.get_N(); }
	plint get_ny() const { return arithmetic_y.get_N(); }
	plint get_nz() const { return arithmetic_z.get_N(); }

	// Remap an int to [0, nx-1], [0, ny-1], [0, nz-1], respectively.
	// If a direction is non-periodic, i is unchanged.
	plint remap_index_x(const plint & i) const
	{
		return arithmetic_x.remap_int(i);
	}

	plint remap_index_x(const plint & i, const plint & offset) const
	{
		return arithmetic_x.remap_int(i, offset);
	}

	plint remap_index_y(const plint & i) const
	{
		return arithmetic_y.remap_int(i);
	}

	plint remap_index_y(const plint & i, const plint & offset) const
	{
		return arithmetic_y.remap_int(i, offset);
	}

	plint remap_index_z(const plint & i) const
	{
		return arithmetic_z.remap_int(i);
	}

	plint remap_index_z(const plint & i, const plint & offset) const
	{
		return arithmetic_z.remap_int(i, offset);
	}

	// Integer difference (positive and in range [0, N])
	plint int_diff_x(plint x1, plint x2) const
	{
		return arithmetic_x.int_diff(x1, x2);
	}

	plint int_diff_y(plint y1, plint y2) const
	{
		return arithmetic_y.int_diff(y1, y2);
	}

	plint int_diff_z(plint z1, plint z2) const
	{
		return arithmetic_z.int_diff(z1, z2);
	}

	// Remap a position to [x0, x0+nx) x [y0, y0+ny) x [z0, z0+nz)
	// A component of the array remains unchanged if the corresponding direction is non-periodic.
	void remap_position(Array<T, 3> & pos) const
	{
		remap_position_x(pos[0]);
		remap_position_y(pos[1]);
		remap_position_z(pos[2]);
	}

	void remap_position_x(T & x) const
	{
		arithmetic_x.remap_pos(x);
	}

	void remap_position_y(T & y) const
	{
		arithmetic_y.remap_pos(y);
	}

	void remap_position_z(T & z) const
	{
		arithmetic_z.remap_pos(z);
	}

	// Push the array pos as close to ref as possible, only adding/subtracting multiples of nx, ny, nz
	// to each component of x, respectively.
	void shift_periodically_to_minimize_distance_to(const Array<T, 3> & ref, Array<T, 3> & pos) const
	{
		arithmetic_x.shift_toward_value(ref[0], pos[0]);
		arithmetic_y.shift_toward_value(ref[1], pos[1]);
		arithmetic_z.shift_toward_value(ref[2], pos[2]);
	}

	// Minimal distance between two values (in modular sense) in a given direction
	T dist_x(const T & p0, const T & p1) const { return arithmetic_x.min_mod_diff(p0, p1); }
	T dist_y(const T & p0, const T & p1) const { return arithmetic_y.min_mod_diff(p0, p1); }
	T dist_z(const T & p0, const T & p1) const { return arithmetic_z.min_mod_diff(p0, p1); }

	// Distance between two 2D vectors
	T dist_sqr_xy(const Array<T, 2> & v1, const Array<T, 2> & v2) const
	{
		return 	detail::pow2(arithmetic_x.min_mod_diff(v1[0], v2[0])) +
				detail::pow2(arithmetic_y.min_mod_diff(v1[1], v2[1]));
	}

	T dist_sqr_xz(const Array<T, 2> & v1, const Array<T, 2> & v2) const
	{
		return 	detail::pow2(arithmetic_x.min_mod_diff(v1[0], v2[0])) +
				detail::pow2(arithmetic_z.min_mod_diff(v1[1], v2[1]));
	}

	T dist_sqr_yz(const Array<T, 2> & v1, const Array<T, 2> & v2) const
	{
		return 	detail::pow2(arithmetic_y.min_mod_diff(v1[0], v2[0])) +
				detail::pow2(arithmetic_z.min_mod_diff(v1[1], v2[1]));
	}

	// Distance between two 3D vectors
	T dist_sqr(const Array<T, 3> & v1, const Array<T, 3> & v2) const
	{
		return 	detail::pow2(arithmetic_x.min_mod_diff(v1[0], v2[0])) +
				detail::pow2(arithmetic_y.min_mod_diff(v1[1], v2[1])) +
				detail::pow2(arithmetic_z.min_mod_diff(v1[2], v2[2]));
	}

	T dist(const Array<T, 3> & v1, const Array<T, 3> & v2) const
	{
		return std::sqrt(dist_sqr(v1, v2));
	}

	Array<T, 3> vec_diff(const Array<T, 3> & v1, const Array<T, 3> & v2) const
	{
		return Array<T, 3>(
				arithmetic_x.min_mod_diff(v1[0], v2[0]),
				arithmetic_y.min_mod_diff(v1[1], v2[1]),
				arithmetic_z.min_mod_diff(v1[2], v2[2])
		);
	}

private:
	Box3D domain;
	detail::Arithmetic1D<T, pX> arithmetic_x;
	detail::Arithmetic1D<T, pY> arithmetic_y;
	detail::Arithmetic1D<T, pZ> arithmetic_z;
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
