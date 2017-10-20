/*
 * TypeDeduction.h
 *
 *  Created on: Sep 16, 2015
 *      Author: niber
 */

#ifndef TYPEDEDUCTION_H_
#define TYPEDEDUCTION_H_

// Methods to strip away pointer or reference attributes from a type.
// e.g. DeduceType<double &>::type is double
template<class T> struct DeduceType { typedef T type; };
template<class T> struct DeduceType<T *> { typedef T type; };
template<class T> struct DeduceType<T &> { typedef T type; };

// Convert a pointer or a reference to a reference
template<class T>
inline T & deref_maybe(T * t)
{
	return *t;
}

template<class T>
inline T & deref_maybe(T & t)
{
	return t;
}



#endif /* TYPEDEDUCTION_H_ */
