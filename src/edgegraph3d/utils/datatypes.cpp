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


#include "datatypes.hpp"

#include <stdexcept>

#include "glm.hpp"

float& vec6::operator[](const int index) {
	switch(index) {
	case 0:
		return x;
	case 1:
		return y;
	case 2:
		return z;
	case 3:
		return i;
	case 4:
		return j;
	case 5:
		return k;
	default:
		std::invalid_argument("Invalid index!");
	}
}

vec6::vec6(const float x,const float y,const float z,const float i,const float j,const float k) : x(x), y(y), z(z), i(i), j(j), k(k) {}

vec6::vec6(const glm::vec3 &s,const glm::vec3 &e) : x(s[0]), y(s[1]), z(s[2]), i(e[0]), j(e[1]), k(e[2]) {}

glm::vec3 vec6::start() {
	return glm::vec3(x,y,z);
}
glm::vec3 vec6::end() {
	return glm::vec3(i,j,k);
}

glm::vec3 vec6::direction() {
	return glm::vec3(i-x,j-y,k-z);
}

bool vec6::operator <(const vec6& v) const
{
    return (x < v.x) || ((!(v.x < x)) && (y < v.y)) || ((!(v.x < x) && !(v.y < y)) && (z < v.z)) || ((!(v.x < x) && !(v.y < y) && !(v.z < z)) && (i < v.i))
    		 || ((!(v.x < x) && !(v.y < y) && !(v.z < z) && !(v.i < i)) && (j < v.j)) || ((!(v.x < x) && !(v.y < y) && !(v.z < z) && !(v.i < i) && !(v.j < j)) && (k < v.k));
}

ray2d::ray2d(const glm::vec2 &start,const glm::vec2 &dir) : start(start), dir(dir) {}

