/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2013 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/** \file
 * Templates for finding indexes for a specific subset of the neighborhood
 *  -- generic implementation.
 */
#ifndef INDEX_UTIL_HH
#define INDEX_UTIL_HH

#include "core/indexUtil.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
Index<T,Descriptor>::Index(plint direction_)
    : direction(direction_)
{
    PLB_ASSERT( direction<Descriptor<T>::d );
}

template<typename T, template<typename U> class Descriptor>
IndexCollection Index<T,Descriptor>::operator==(plint value) const {
    std::vector<plint> indexes;
    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        if (Descriptor<T>::c[iPop][direction]==value) {
            indexes.push_back(iPop);
        }
    }
    return IndexCollection(indexes);
}

template<typename T, template<typename U> class Descriptor>
IndexCollection Index<T,Descriptor>::operator<(plint value) const {
    std::vector<plint> indexes;
    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        if (Descriptor<T>::c[iPop][direction]<value) {
            indexes.push_back(iPop);
        }
    }
    return IndexCollection(indexes);
}

template<typename T, template<typename U> class Descriptor>
IndexCollection Index<T,Descriptor>::operator<=(plint value) const {
    std::vector<plint> indexes;
    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        if (Descriptor<T>::c[iPop][direction]<=value) {
            indexes.push_back(iPop);
        }
    }
    return IndexCollection(indexes);
}

template<typename T, template<typename U> class Descriptor>
IndexCollection Index<T,Descriptor>::operator>(plint value) const {
    std::vector<plint> indexes;
    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        if (Descriptor<T>::c[iPop][direction]>value) {
            indexes.push_back(iPop);
        }
    }
    return IndexCollection(indexes);
}

template<typename T, template<typename U> class Descriptor>
IndexCollection Index<T,Descriptor>::operator>=(plint value) const {
    std::vector<plint> indexes;
    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        if (Descriptor<T>::c[iPop][direction]>=value) {
            indexes.push_back(iPop);
        }
    }
    return IndexCollection(indexes);
}

}  // namespace plb

#endif  // INDEX_UTIL_HH