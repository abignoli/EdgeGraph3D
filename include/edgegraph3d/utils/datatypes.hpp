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


#ifndef INCLUDE_EDGEGRAPH3D_UTILS_DATATYPES_HPP_
#define INCLUDE_EDGEGRAPH3D_UTILS_DATATYPES_HPP_

#include <set>

#include "segment_edge_manager.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"

/**
 * correspondences -> vector of starting image is empty (i.e. starting intersection not included)
 */
typedef struct potential_match {
	SegmentEdgeManager::segment_point starting_intersection;
	std::vector<std::vector<uint>> correspondences;
	vector<SegmentEdgeManager::segment_point> correspondences_data;
} potential_match;

typedef struct vec6 {
	float x;
	float y;
	float z;
	float i;
	float j;
	float k;
	float& operator[](const int index);
	vec6::vec6(const float x,const float y,const float z,const float i,const float j,const float k);
	vec6(const glm::vec3 &s,const glm::vec3 &e);
	glm::vec3 start();
	glm::vec3 end();
	glm::vec3 direction();

    bool operator <(const vec6& v) const;
} vec6;

typedef struct ray2d {
	glm::vec2 start;
	glm::vec2 dir;
	ray2d(const glm::vec2 &start,const glm::vec2 &dir);
} ray2d;

struct lex_compare {
    bool operator() (const glm::vec3& a, const glm::vec3& b) const{
    	return a.x == b.x && a.y == b.y && a.z == b.z;
    }
};

typedef set<glm::vec3, lex_compare> set3dpoints;


#endif /* INCLUDE_EDGEGRAPH3D_UTILS_DATATYPES_HPP_ */
