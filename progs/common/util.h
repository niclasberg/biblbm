/*
 * util.h
 *
 *  Created on: Jul 21, 2015
 *      Author: niber
 */

#ifndef UTIL_H_
#define UTIL_H_
#include "core/globalDefs.h"
#include "parallelism/mpiManager.h"
#include <cstdlib>

template<class T, plb::pluint N>
void random_vector(plb::Array<T, N> & res)
{
	for(plb::pluint i = 0; i < N; ++i)
		res[i] = std::rand() / RAND_MAX;
}


#endif /* UTIL_H_ */
