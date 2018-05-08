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


#ifndef INCLUDE_EDGEGRAPH3D_UTILS_GEOMETRY_TRIANGULATION_HPP_
#define INCLUDE_EDGEGRAPH3D_UTILS_GEOMETRY_TRIANGULATION_HPP_

#include <vector>

#include "polyline_graph_2d.hpp"
#include "polyline_graph_2d_hmap_impl.hpp"
#include "polyline_2d_map_search.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"
#include "global_switches.hpp"
#include "glm.hpp"

using namespace std;

#define MAX_3DPOINT_PROJECTIONDISTSQ_EXPANDALLVIEWS 16.0

void compute_3d_point_coords(const SfMData &sfmd,
		const vector<glm::vec2> &selected_2d_reprojections_coords,
		const vector<int> &selected_2d_reprojections_ids,
		glm::vec3 &new_point_data,
		bool &valid);

void compute_all_potential_3d_points_3views(const SfMData &sfmd,
		const tuple<vector<PolyLineGraph2D::plg_point>,vector<PolyLineGraph2D::plg_point>,vector<PolyLineGraph2D::plg_point>> &all_2d_reprojections_3views,
		const tuple<int,int,int> &all_2d_reprojections_ids_3views,
		vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &potential_new_points_data);

void compute_all_potential_3d_points_3views_plg_following_newpoint_compatibility(const SfMData &sfmd,const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs,
		const tuple<vector<PolyLineGraph2D::plg_point>,vector<PolyLineGraph2D::plg_point>,vector<PolyLineGraph2D::plg_point>> &all_2d_reprojections_3views,
		const tuple<int,int,int> &all_2d_reprojections_ids_3views,
		vector<new_3dpoint_and_sides_plgp_matches> &potential_new_points_and_following_data_directions);

void compute_3D_point_multiple_views_plg_following(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const int starting_plg_id, const vector<vector<PolyLineGraph2D::plg_point>> &epipolar_correspondences, new_3dpoint_and_sides_plgp_matches &p3d_with_sides, bool &valid);

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> compute_3D_point_multiple_views_plg_following_vecpoints(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const int starting_plg_id, const vector<vector<PolyLineGraph2D::plg_point>> &epipolar_correspondences);

void compute_3d_point_coords_combinations(const SfMData &sfmd,
		const vector<glm::vec2> &all_2d_reprojections_coords,
		const vector<int> &all_2d_reprojections_ids,
		const int min_combinations,
		vector<glm::vec2> &selected_2d_reprojections_coords,
		vector<int> &selected_2d_reprojections_ids,
		vector<bool> &selected,
		glm::vec3 &new_point_data,
		bool &valid);

void em_estimate3Dpositions(const SfMData &sfmd, const vector<PolyLineGraph2D::plg_point> &selected_2d_reprojections_coords, const vector<int> &selected_2d_reprojections_ids, glm::vec3 &triangulated_point, bool &valid);

void em_estimate3Dpositions(const SfMData &sfmd, const vector<glm::vec2> &selected_2d_reprojections_coords, const vector<int> &selected_2d_reprojections_ids, glm::vec3 &triangulated_point, bool &valid);

bool compatible_new_observation_to_3Dpositions(const SfMData &sfmd, const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &current_point, const glm::vec2 new_coords, const int new_viewpoint_id, glm::vec3 &triangulated_point);

bool compatible_new_observation_to_3Dpositions(const SfMData &sfmd, const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &current_point, const PolyLineGraph2D::plg_point &new_plgp, const int new_viewpoint_id, glm::vec3 &triangulated_point);
bool compatible_new_observation_to_3Dpositions_update(const SfMData &sfmd, std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &current_point, const PolyLineGraph2D::plg_point &new_plgp, const int new_viewpoint_id);

void em_add_new_observation_to_3Dpositions(const SfMData &sfmd, const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &current_point, const glm::vec2 new_coords, const int new_viewpoint_id, glm::vec3 &triangulated_point, bool &valid);
void em_add_new_observation_to_3Dpositions(const SfMData &sfmd, const glm::vec3 &current_point_coords, const vector<glm::vec2> &current_point_observation_coords, const vector<int> &current_point_observation_ids, const glm::vec2 new_coords, const int new_viewpoint_id, glm::vec3 &triangulated_point, bool &valid);

void compute_3d_point(const SfMData &sfmd,
		const vector<glm::vec2> &selected_2d_reprojections_coords,
		const vector<int> &selected_2d_reprojections_ids,
		std::tuple<glm::vec3, vector<glm::vec2>, vector<int>> &new_point_data,
		bool &valid);

void expand_point_to_other_views_expandallviews(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const int starting_plg_id, const vector<vector<PolyLineGraph2D::plg_point>> &epipolar_correspondences, const int selected_views[], const vector<PolyLine2DMapSearch> &plmaps, new_3dpoint_and_sides_plgp_matches &current_p3d_with_sides);
void compute_3D_point_multiple_views_plg_following_expandallviews(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const int starting_plg_id, const vector<vector<PolyLineGraph2D::plg_point>> &epipolar_correspondences, const vector<PolyLine2DMapSearch> &plmaps, new_3dpoint_and_sides_plgp_matches &p3d_with_sides, bool &valid);
vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> compute_3D_point_multiple_views_plg_following_expandallviews_vector(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const int starting_plg_id, const vector<vector<PolyLineGraph2D::plg_point>> &epipolar_correspondences, const vector<PolyLine2DMapSearch> &plmaps);
vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> compute_3D_point_multiple_views_plg_following_vecpoints_expandallviews(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const int starting_plg_id, const vector<vector<PolyLineGraph2D::plg_point>> &epipolar_correspondences, const vector<PolyLine2DMapSearch> &plmaps);

#endif /* INCLUDE_EDGEGRAPH3D_UTILS_GEOMETRY_TRIANGULATION_HPP_ */
