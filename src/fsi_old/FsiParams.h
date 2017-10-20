/*
 * FsiParams.h
 *
 *  Created on: Apr 17, 2014
 *      Author: niber
 */

#ifndef FSIPARAMS_H_
#define FSIPARAMS_H_
#include

namespace plb {

namespace fsi {

template<class T>
class FsiParams {
public:
	FsiParams();

private:
	//
	T potential_strength_nondim, potential_strength_lu;


	IncomprFlowParam<T> flow_params;
};

}

}

#endif /* FSIPARAMS_H_ */
