/*
 ***********************************************************************
 *
 * 							  EdgeGraph3D
 *
 *         Copyright (C) 2018 Andrea Bignoli (andrea.bignoli@gmail.com)
 *                         All rights reserved
 *
 ***********************************************************************
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 *
*/


#ifndef INCLUDE_EDGEGRAPH3D_UTILS_SERIALIZATION_UTILITIES_HPP_
#define INCLUDE_EDGEGRAPH3D_UTILS_SERIALIZATION_UTILITIES_HPP_


#include "glm.hpp"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>

#include "global_defines.hpp"
#include "global_switches.hpp"

namespace boost {
namespace serialization {

template<class Archive>
void serialize(Archive & ar, glm::vec2 & v, const unsigned int version)
{
    ar & v.x;
    ar & v.y;
}

template<class Archive>
void serialize(Archive & ar, glm::vec3 & v, const unsigned int version)
{
    ar & v.x;
    ar & v.y;
    ar & v.z;
}

} // namespace serialization
} // namespace boost


#endif /* INCLUDE_EDGEGRAPH3D_UTILS_SERIALIZATION_UTILITIES_HPP_ */
