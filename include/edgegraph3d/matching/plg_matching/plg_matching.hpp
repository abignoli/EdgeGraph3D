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


#ifndef INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_PLG_MATCHING_HPP_
#define INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_PLG_MATCHING_HPP_

#include <utility>

#include "polyline_graph_2d.hpp"
#include "glm.hpp"

#define PLG_FOLLOW_FIRST_IMAGE_DISTANCE 10.0
#define PLG_FOLLOW_CORRESPONDENCE_IMAGE_DISTANCE_MIN (PLG_FOLLOW_FIRST_IMAGE_DISTANCE/2)
#define PLG_FOLLOW_CORRESPONDENCE_IMAGE_DISTANCE_MAX (PLG_FOLLOW_FIRST_IMAGE_DISTANCE*2)

#include "datatypes.hpp"
#include "triangulation.hpp"
#include <iostream>
#include <vector>
#include <tuple>
#include "edge_graph_3d_utilities.hpp"
#include "geometric_utilities.hpp"
#include "drawing_utilities.hpp"
#include <string>
#include "convert_edge_images_pixel_to_segment.hpp"
#include "test_utilities.hpp"
#include "polyline_graph_2d_hmap_impl.hpp"
#include "plg_edge_manager_closest_only.hpp"

#include "global_switches.hpp"
#include "global_defines.hpp"

#define SWITCH_PLG_MATCHING_ADDPOINT_BOTHDIR_ONE

#define PLG_MATCHING_TRIANGULATION_MINIMUM_AMOUNT_OF_POINTS 3

vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> follow_plgs_from_match(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs, const std::tuple<glm::vec3,vector<PolyLineGraph2D::plg_point>,vector<int>> &matches, bool &valid);

pair<vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>>,vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>>> follow_plgs_from_match2(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs, const std::tuple<glm::vec3,vector<PolyLineGraph2D::plg_point>,vector<int>> &matches, bool &valid);

pair<vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>>,vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>>> follow_plgs_from_match3(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs, const std::tuple<glm::vec3,vector<PolyLineGraph2D::plg_point>,vector<int>> &matches, bool &valid);

void follow_plgs_from_match4(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs, const new_3dpoint_plgp_matches &matches, vector<ulong> &directions1, bool &direction1_valid, vector<new_3dpoint_plgp_matches> &valid_points_direction1, vector<ulong> &directions2, bool &direction2_valid, vector<new_3dpoint_plgp_matches> &valid_points_direction2);
bool compatible_new_plg_point(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs, const new_3dpoint_plgp_matches &matches, vector<ulong> &directions1, bool &direction1_valid, vector<new_3dpoint_plgp_matches> &valid_points_direction1, vector<ulong> &directions2, bool &direction2_valid, vector<new_3dpoint_plgp_matches> &valid_points_direction2);
bool add_view_to_3dpoint_and_sides_plgp_matches(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs, new_3dpoint_and_sides_plgp_matches &cur_pts, const int current_plg_id, const PolyLineGraph2D::plg_point &current_plgp);
pair<int,int> add_view_to_3dpoint_and_sides_plgp_matches_vector(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs, vector<new_3dpoint_plgp_matches> &cur_pts, vector<ulong> &start_dirs, vector<ulong> &end_dirs, const int current_plg_id, const PolyLineGraph2D::plg_point &current_plgp, int start_interval_index, int cur_point_index, int end_interval_index, bool &success);
#endif /* INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_PLG_MATCHING_HPP_ */
