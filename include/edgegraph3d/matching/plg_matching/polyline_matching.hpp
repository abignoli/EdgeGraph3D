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


#ifndef INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_POLYLINE_MATCHING_HPP_
#define INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_POLYLINE_MATCHING_HPP_


#include <set>
#include <vector>

#include "polyline_graph_2d_hmap_impl.hpp"
#include "polyline_2d_map_search.hpp"
#include "polyline_graph_2d.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"
#include "global_switches.hpp"
#include "glm.hpp"

class PLGMatchesManager;

using namespace std;

#define SPLIT_INTERVAL_DISTANCE 20.0

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> find_new_3d_points_from_compatible_polylines(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const vector<set<ulong>> &potentially_compatible_polylines, PLGMatchesManager &plgmm);
//vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> find_new_3d_points_from_compatible_polylines_parallel(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const vector<set<ulong>> &potentially_compatible_polylines, PLGMatchesManager &plgmm);
vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> find_new_3d_points_from_compatible_polylines_expandallviews(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const vector<set<ulong>> &potentially_compatible_polylines, PLGMatchesManager &plgmm, const vector<PolyLine2DMapSearch> &plmaps);
vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> find_new_3d_points_from_compatible_polylines_expandallviews_parallel(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const vector<set<ulong>> &potentially_compatible_polylines, PLGMatchesManager &plgmm, const vector<PolyLine2DMapSearch> &plmaps);

#endif /* INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_POLYLINE_MATCHING_HPP_ */
