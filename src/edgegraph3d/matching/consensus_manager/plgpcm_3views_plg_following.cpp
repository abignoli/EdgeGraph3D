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

#include "SfMData.h"
#include "plgpcm_3views_plg_following.hpp"

#include "glm.hpp"

#include "triangulation.hpp"
#include "global_switches.hpp"

PLGPCM3ViewsPLGFollowing::PLGPCM3ViewsPLGFollowing(vector<Mat> &input_imgs, const SfMData &input_sfmd, const Size &input_sz, const Mat** all_fundamental_matrices, vector<PolyLineGraph2DHMapImpl> &plgs, const vector<PolyLine2DMapSearch> &plmaps) : PLGPConsensusManager(input_imgs, input_sfmd, input_sz, all_fundamental_matrices, plgs), plmaps(plmaps) {}

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> PLGPCM3ViewsPLGFollowing::consensus_strategy_single_point_single_intersection(const int starting_img_id, const int refpoint, PolyLineGraph2D::plg_point intersection_on_starting_img,std::vector<std::vector<PolyLineGraph2D::plg_point>> intersection_correspondences) {
	std::vector<std::vector<PolyLineGraph2D::plg_point>> intersection_correspondences_all(plgs.size());
	for(int i=0; i < sfmd.camViewingPointN_[refpoint].size(); i++)
		intersection_correspondences_all[sfmd.camViewingPointN_[refpoint][i]] = intersection_correspondences[i];

#if defined(GLOBAL_SWITCH_USEEXPANDALLVIEWSVECTOR)
	return compute_3D_point_multiple_views_plg_following_expandallviews_vector(sfmd, plgs, all_fundamental_matrices, starting_img_id, intersection_correspondences_all,plmaps);
#else
	return compute_3D_point_multiple_views_plg_following_vecpoints_expandallviews(sfmd, plgs, all_fundamental_matrices, starting_img_id, intersection_correspondences_all,plmaps);
#endif
}

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> PLGPCM3ViewsPLGFollowing::consensus_strategy_single_point(const int starting_img_id, const int refpoint, const std::pair<std::vector<PolyLineGraph2D::plg_point>,std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>>> &intersections_and_correspondences_pair) {
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> res;
	for(int i=0; i < intersections_and_correspondences_pair.first.size(); i++) {
		vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> cur_res = consensus_strategy_single_point_single_intersection(starting_img_id, refpoint, intersections_and_correspondences_pair.first[i],intersections_and_correspondences_pair.second[i]);
		for(auto &newp3d : cur_res)
			res.push_back(newp3d);
	}
	return res;
}

vector<vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>>> PLGPCM3ViewsPLGFollowing::consensus_strategy_single_point_vector(const int starting_img_id, const int refpoint, const std::pair<std::vector<PolyLineGraph2D::plg_point>,std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>>> &intersections_and_correspondences_pair) {
	vector<vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>>> res;
	for(int i=0; i < intersections_and_correspondences_pair.first.size(); i++) {
		vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> cur_res = consensus_strategy_single_point_single_intersection(starting_img_id, refpoint, intersections_and_correspondences_pair.first[i],intersections_and_correspondences_pair.second[i]);
		res.push_back(cur_res);
	}
	return res;
}
