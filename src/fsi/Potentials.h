#ifndef POTENTIALS_H_
#define POTENTIALS_H_

namespace plb {

namespace fsi {

template<class T>
class LennardJones12_6 {
public:
	LennardJones12_6(T r_min_, T A_) : r_min(r_min_), A(A_), c() { }

	T operator()(T dist) const
	{
		const T x = std::pow(r_min/dist, 6);
		//std::cout << 12.0*A/dist * (x*x - x) << std::endl;
		return 12.0*A/dist * (x*x - x);
	}

	T get_repulsion_distance() const { return r_min; }

private:
	T c;
	T r_min;
	T A;
};

template<class T>
class MorsePotential {
public:
	MorsePotential(T r_min_, T beta_, T strength_)
	: r_min(r_min_), beta(beta_), strength(strength_), two_beta_strenght(2*beta*strength)
	{
		r_cut = 2*r_min;
	}

	T operator()(T dist) const
	{
		return eval(dist);
	}

	T get_cutoff_distance() const { return r_min; }
	T get_repulsion_distance() const { return r_min; }

private:
	T eval(T dist) const
	{
		if(dist < 0)
			dist = (T) 0;
		const T x = std::exp(beta*(r_min-dist));
		return two_beta_strenght * x*(x - 1);
	}

	T r_min, beta, strength, two_beta_strenght, r_cut;
};

template<class T>
class SpringPotential {
public:
	SpringPotential(T min_dist_, T strength_) : min_dist(min_dist_), strength(strength_) { }
	T operator()(T dist) const
	{
		return (min_dist - dist) * strength;
	}

	T get_repulsion_distance() const { return min_dist; }
	T get_cutoff_distance() const { return min_dist; }

private:
	T min_dist, strength;
};

}

}




#endif /* POTENTIALS_H_ */
