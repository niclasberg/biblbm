/*
 * utils.h
 *
 *  Created on: Jun 12, 2015
 *      Author: niber
 */

#ifndef UTILS_H_
#define UTILS_H_
#include <vector>
#include <algorithm>
#include <ostream>
#include <istream>

namespace plb {

namespace fsi {

namespace utils {

template<class T>
void pack(std::vector<char> & vec, const T & data)
{
	const char * const buf = reinterpret_cast<const char *>(&data);
	vec.insert(vec.end(), buf, buf + sizeof(T));
}

template<class T>
void pack(char * vec, const T & data)
{
	const char * const buf = reinterpret_cast<const char *>(&data);
	std::copy(buf, buf + sizeof(T), vec);
}

template<class T>
void pack(std::ostream & buffer, const T & data)
{
	buffer.write(reinterpret_cast<const char *>(&data), sizeof(T));
}

template<class T>
void unpack(std::vector<char>::iterator & it, T & data)
{
	std::copy(it, it+sizeof(T), reinterpret_cast<char *>(&data));
	it += sizeof(T);
}

template<class T>
void unpack(char ** it, T & data)
{
	std::copy(*it, (*it)+sizeof(T), reinterpret_cast<char *>(&data));
	*it += sizeof(T);
}

template<class T>
void unpack(char *& it, T & data)
{
	std::copy(it, it+sizeof(T), reinterpret_cast<char *>(&data));
	it += sizeof(T);
}

template<class T>
void unpack(std::istream & buffer, T & data)
{
	char * d = reinterpret_cast<char *>(&data);
	buffer.read(d, sizeof(T));
}

}

std::ostream & operator<<(std::ostream & s, const Box3D & box)
{
	s << "( [" << box.x0 << ", " << box.x1 << "], [" << box.y0 << ", " << box.y1 << "], [" << box.z0 << ", " << box.z1 << "] )";
	return s;
}

}

}

#endif /* UTILS_H_ */
