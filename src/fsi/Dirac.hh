#ifndef DIRAC_HH_
#define DIRAC_HH_
#include "Dirac.h"

namespace plb {

namespace fsi {

namespace detail {

// Optimize so that the weights are as close to 1 as possible, and constrain the optimization
// so that the final Dirac function fulfills the moment conditions (sum dirac(p) = 1, sum x*dirac(p) = sum y*dirac(p) = sum z*dirac(p) = 0)
double dirac_weights_cost_function(const std::vector<double> & x, std::vector<double> & grad, void * data)
{
	if(!grad.empty()) {
		for(int i=0; i < grad.size(); ++i)
			grad[i] = 2.*(x[i] - 1);
	}
	double ret = 0.;
	for(int i=0; i < x.size(); ++i)
		ret += (x[i] - 1.)*(x[i] - 1.);
	return ret;
}

template<class Dirac>
double constrain_zero_moment(const std::vector<double> & x, std::vector<double> & grad, void * data)
{
	std::vector<typename DiracWithMissingPoints<double, Dirac>::DiracPoint> & dirac_points =
			*reinterpret_cast<std::vector<typename DiracWithMissingPoints<double, Dirac>::DiracPoint> *>(data);
	double ret = -1.;
	for(plint i = 0; i < x.size(); ++i)
		ret += dirac_points[i].dirac_val * x[i];
	if(!grad.empty()) {
		for(plint i=0; i < grad.size(); ++i)
			grad[i] = dirac_points[i].dirac_val;
	}
	return ret;
}

template<class Dirac, int dim>
double constrain_first_moment(const std::vector<double> & x, std::vector<double> & grad, void * data)
{
	std::vector<typename DiracWithMissingPoints<double, Dirac>::DiracPoint> & dirac_points =
			*reinterpret_cast<std::vector<typename DiracWithMissingPoints<double, Dirac>::DiracPoint> *>(data);
	double ret = 0.;
	for(plint i = 0; i < x.size(); ++i)
		ret += dirac_points[i].dirac_val * x[i] * dirac_points[i].dx[dim];
	if(!grad.empty()) {
		for(plint i=0; i < grad.size(); ++i)
			grad[i] = dirac_points[i].dirac_val * dirac_points[i].dx[dim];
	}
	return ret;
}

} /* namespace detail */

template<class T, class Dirac>
void DiracWithMissingPoints<T, Dirac>::computeWeights() {
	// Create nodes
	_dirac_points.clear();
	plint huge = std::numeric_limits<plint>::max();
	Box3D bb(huge, -huge, huge, -huge, huge, -huge);
	for(plint i = 0; i < Dirac::support; ++i)
		for(plint j = 0; j < Dirac::support; ++j)
			for(plint k = 0; k < Dirac::support; ++k)
				if(_nodeIsValid[i][j][k]) {
					DiracPoint p;
					p.node_pos = Dot3D(i+_i0.x, j+_i0.y, k+_i0.z);
					p.dx = _x0 - Array<T, 3>((T)p.node_pos.x, (T)p.node_pos.y, (T)p.node_pos.z);
					p.dirac_val = Dirac::eval(p.dx[0], p.dx[1], p.dx[2]);
					p.weight = 1.;
					_dirac_points.push_back(p);
					bb.x0 = std::min(bb.x0, i);
					bb.x1 = std::max(bb.x1, i); 
					bb.y0 = std::min(bb.y0, j);
					bb.y1 = std::max(bb.y1, j);
					bb.z0 = std::min(bb.z0, k);
					bb.z1 = std::max(bb.z1, k);
				}

	// If no points are missing, we can just keep the weights at 1 and use the ordinary Dirac function for
	// interpolation
	if(bb.getNx() == Dirac::support && bb.getNy() == Dirac::support && bb.getNz() == Dirac::support)
		return;

	if(_dirac_points.empty()) {
		// No points were found
		std::cerr << "Warning: trying to interpolate to a particle node outside of the domain" << std::endl
				<< "Node position: (" << _x0[0] << ", " << _x0[1] << ", " << _x0[2] << ")" << std::endl;
	} else {
		// Setup optimization problem
		// It is not possible to fulfill the moment condition in a direction if
		// less than two points are included in that direction.
		/*nlopt::opt opt(nlopt::LD_SLSQP, _dirac_points.size());
		opt.set_min_objective(detail::dirac_weights_cost_function, 0);
		opt.add_equality_constraint(detail::constrain_zero_moment<Dirac>, reinterpret_cast<void *>(&_dirac_points), 1e-6);
		if(bb.x1 - bb.x0 > 0)
			opt.add_equality_constraint(detail::constrain_first_moment<Dirac, 0>, reinterpret_cast<void *>(&_dirac_points), 1e-6);
		if(bb.y1 - bb.y0 > 0)
			opt.add_equality_constraint(detail::constrain_first_moment<Dirac, 1>, reinterpret_cast<void *>(&_dirac_points), 1e-6);
		if(bb.z1 - bb.z0 > 0)
			opt.add_equality_constraint(detail::constrain_first_moment<Dirac, 2>, reinterpret_cast<void *>(&_dirac_points), 1e-6);
		opt.set_xtol_rel(1e-4);

		// Start guess: all weights are 1
		std::vector<T> weights(_dirac_points.size(), 1.);
		T minf;
		//try {
			nlopt::result result = opt.optimize(weights, minf);
		//} catch()

		for(int i = 0; i < weights.size(); ++i)
			_dirac_points[i].weight = weights[i];*/

		T dirac_sum = 0.;
		for(int i = 0; i < _dirac_points.size(); ++i)
			dirac_sum += _dirac_points[i].dirac_val;
		for(int i = 0; i < _dirac_points.size(); ++i)
			_dirac_points[i].weight = 1. / dirac_sum;
	}
}

} /* namespace fsi */

}/* namespace plb */


#endif /* DIRAC_HH_ */
