#ifndef DIRAC_H_
#define DIRAC_H_
#include "core/array.h"

namespace plb {

namespace fsi {

template<class T, int N>
class PeskinDirac { };

template<class T>
struct PeskinDirac<T, 3> {
	static T half_support;
	static pluint support;

	static T eval(T x)
	{
		T dx = std::fabs(x);
		return 0.125 * (dx <= 1 ? (3 - 2*dx + std::sqrt(1 + 4*dx*(1-dx))) : (5 - 2*dx - std::sqrt(-7 + 4*dx*(3 - dx))));
	}

	static T eval(T x, T y, T z)
	{
		T dx = std::fabs(x), dy = std::fabs(y), dz = std::fabs(z);
		return 0.001953125 *
				(dx <= 1 ? (3 - 2*dx + std::sqrt(1 + 4*dx*(1-dx))) : (5 - 2*dx - std::sqrt(-7 + 4*dx*(3 - dx)))) *
				(dy <= 1 ? (3 - 2*dy + std::sqrt(1 + 4*dy*(1-dy))) : (5 - 2*dy - std::sqrt(-7 + 4*dy*(3 - dy)))) *
				(dz <= 1 ? (3 - 2*dz + std::sqrt(1 + 4*dz*(1-dz))) : (5 - 2*dz - std::sqrt(-7 + 4*dz*(3 - dy))));
	}
};

template<class T> pluint PeskinDirac<T, 3>::support = 4;
template<class T> T PeskinDirac<T, 3>::half_support = 2;

template<class T, int N>
class RomaDirac { };

template<class T>
struct RomaDirac<T, 3> {
	static T half_support;
	static pluint support;

	static T eval(T x)
	{
		T dx = std::fabs(x);
		PLB_PRECONDITION(dx <= 1.5)
		return (1.0 / 6.0) * (dx > 0.5 ? (5.0 - 3*dx - std::sqrt(-3*(1-dx)*(1-dx)+1)) : 2.0 * (1 + std::sqrt(1-3*dx*dx)));
	}

	static T eval(T x, T y, T z)
	{
		// For performance: assume that x, y, z is in the interval [-1.5, 1.5]
		T dx = std::fabs(x), dy = std::fabs(y), dz = std::fabs(z);
		return (1.0 / 216.0) *
				(dx > 0.5 ? (5.0 - 3*dx - std::sqrt(-3*(1-dx)*(1-dx)+1)) : 2.0 * (1 + std::sqrt(1-3*dx*dx))) *
				(dy > 0.5 ? (5.0 - 3*dy - std::sqrt(-3*(1-dy)*(1-dy)+1)) : 2.0 * (1 + std::sqrt(1-3*dy*dy))) *
				(dz > 0.5 ? (5.0 - 3*dz - std::sqrt(-3*(1-dz)*(1-dz)+1)) : 2.0 * (1 + std::sqrt(1-3*dz*dz)));
	}
};

template<class T> pluint RomaDirac<T, 3>::support = 3;
template<class T> T RomaDirac<T, 3>::half_support = 1.5;

template<class T, class Dirac, int N>
class SampledDirac { };

template<class T, class Dirac>
struct SampledDirac<T, Dirac, 3> {
public:
	SampledDirac() : values(1024), dx(Dirac::support / (T) 1023), dx_inv((T) 1023 / Dirac::support)
	{
		for(int i = 0; i < 1024; ++i) {
			values[i] = Dirac::eval(-Dirac::half_support + dx*i);
		}
	}

	T eval(T x)
	{
		return values[std::floor(0.5 + (x+Dirac::half_support)*dx_inv)];
	}

	T eval(T x, T y, T z)
	{
		return eval(x)*eval(y)*eval(z);
	}

private:
	std::vector<T> values;
	T dx;
	T dx_inv;
};

template<typename T, typename Dirac>
void get_dirac_compact_support_box(const Array<T, 3> & pos, Box3D & ret)
{
	ret.x0 = std::ceil(pos[0] - Dirac::half_support);
	ret.x1 = ret.x0 + Dirac::support - 1;
	ret.y0 = std::ceil(pos[1] - Dirac::half_support);
	ret.y1 = ret.y0 + Dirac::support - 1;
	ret.z0 = std::ceil(pos[2] - Dirac::half_support);
	ret.z1 = ret.z0 + Dirac::support - 1;
}

} /* namespace fsi */

} /* namespace plb */

#endif /* DIRAC_H_ */
