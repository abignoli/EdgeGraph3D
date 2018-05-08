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


#ifndef INCLUDE_EDGEGRAPH3D_UTILS_PLG_ANALYSIS_UTILITIES_HPP_
#define INCLUDE_EDGEGRAPH3D_UTILS_PLG_ANALYSIS_UTILITIES_HPP_

#include <vector>

#include "polyline_graph_2d_hmap_impl.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"
#include "global_switches.hpp"
#include "glm.hpp"

#define MAX_EDGEREFPOINT_DIST 1.0
#define MAX_EDGEREFPOINT_DISTSQ (MAX_EDGEREFPOINT_DIST*MAX_EDGEREFPOINT_DIST)

int is_edgerefpoint(const vector<vector<glm::vec4>> &vec_segments, const SfMData &sfmd, int refpoint_id);

std::vector<int> find_edgerefpoints(const vector<PolyLineGraph2DHMapImpl> &plgs, const SfMData &sfmd);

#endif /* INCLUDE_EDGEGRAPH3D_UTILS_PLG_ANALYSIS_UTILITIES_HPP_ */
