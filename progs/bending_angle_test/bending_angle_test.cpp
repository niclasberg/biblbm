#include "fsi/TriangleUtils.h"
#include <iostream>

using namespace plb;

int main(int argc, char * argv[])
{
	// Triangle vertices
	Array<Array<double, 3>, 4> v;
	v[0] = Array<double, 3>(0, -1, 0);
	v[1] = Array<double, 3>(1, 0, 0);
	v[2] = Array<double, 3>(0, 1, 0);
	v[3] = Array<double, 3>(-1, 0, 0);
	Array<Array<double, 3>, 4> v2;

	double cos_angle, sin_angle;
	double cos_angle2, sin_angle2;

	Array<Array<double, 3>, 4> g;
	double dtheta = 0.05;
	double dx = 1e-5;

	std::cout << std::scientific;

	for(int itheta = -3; itheta <= 3; ++itheta) {
		double theta = itheta * dtheta;
		v[1][0] = std::cos(theta);
		v[1][2] = -std::sin(theta);
		fsi::tri::grad_angle_between_pair(v[0], v[1], v[2], v[3],
				cos_angle, sin_angle,
				g[0], g[1], g[2], g[3]);
		std::cout << "theta = " << theta << ", cos(angle) = " << cos_angle << " (" << std::cos(theta)
				  << "), sin(angle) = " << sin_angle << " (" << std::sin(theta) << ")" << std::endl;
		/*if(theta != 0.) {
			double denom = 1. + (v1[2] / v1[0])*(v1[2] / v1[0]);

			std::cout << "  grad_x1 = " << g1[0] << " (" << v1[2]/(v1[0]*v1[0]) / denom << "), "
					  << g1[1] << " (0), "
					  << g1[2] << " (" << -1. / v1[0] / denom << ")" << std::endl;

		}*/
		for(int i = 0; i < 4; ++i) {
			std::cout << "  ";
			for(int dim = 0; dim < 3; ++dim) {
				// Perturb the vertex
				v2 = v;
				v2[i][dim] += dx;

				fsi::tri::cos_sin_of_angle_between_pair(v2[0], v2[1], v2[2], v2[3],
						cos_angle2, sin_angle2);
				std::cout << g[i][dim] << " (" << (std::asin(sin_angle2) - std::asin(sin_angle)) / dx << "), ";
			}
			std::cout << std::endl;
		}
	}
}

