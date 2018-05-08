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


#ifndef INCLUDE_EDGEGRAPH3D_FILTERING_FILTERING_CLOSE_PLGPS_HPP_
#define INCLUDE_EDGEGRAPH3D_FILTERING_FILTERING_CLOSE_PLGPS_HPP_

#include <opencv2/core/types.hpp>
#include <vector>

#include "polyline_graph_2d.hpp"
#include "polyline_graph_2d_hmap_impl.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"
#include "glm.hpp"


vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> filter_3d_points_close_2d_array(const vector<PolyLineGraph2DHMapImpl> &plgs, const Size &imgsz, const vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &p3ds);

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> filter3dpoints_close2d(const vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &p3ds);

#endif /* INCLUDE_EDGEGRAPH3D_FILTERING_FILTERING_CLOSE_PLGPS_HPP_ */
