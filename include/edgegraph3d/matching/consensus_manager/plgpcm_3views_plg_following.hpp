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


#ifndef INCLUDE_ALL_COMBINATIONS_CONSENSUS_MANAGER_PLGPCM3VIEWSPLGFOLLOWING_HPP_
#define INCLUDE_ALL_COMBINATIONS_CONSENSUS_MANAGER_PLGPCM3VIEWSPLGFOLLOWING_HPP_

#include <opencv2/core/mat.hpp>
#include <opencv2/core/types.hpp>
#include <utility>
#include <vector>

#include "polyline_graph_2d_hmap_impl.hpp"
#include "plgp_consensus_manager.hpp"
#include "polyline_2d_map_search.hpp"
#include "polyline_graph_2d.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"
#include "global_switches.hpp"
#include "glm.hpp"

/*
 * Manages consensus problem to enable inlier PLGPs detection during simultaneous 3-View PLG exploration
 */
class PLGPCM3ViewsPLGFollowing : public PLGPConsensusManager {
private:
	const vector<PolyLine2DMapSearch> &plmaps;
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> consensus_strategy_single_point_single_intersection(const int starting_img_id, const int refpoint, PolyLineGraph2D::plg_point intersection_on_starting_img,std::vector<std::vector<PolyLineGraph2D::plg_point>> intersection_correspondences);
	vector<vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>>> consensus_strategy_single_point_vector(const int starting_img_id, const int refpoint, const std::pair<std::vector<PolyLineGraph2D::plg_point>,std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>>> &intersections_and_correspondences_pair);
public:
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> consensus_strategy_single_point(const int starting_img_id, const int refpoint, const std::pair<std::vector<PolyLineGraph2D::plg_point>,std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>>> &intersections_and_correspondences_pair);
	PLGPCM3ViewsPLGFollowing(vector<Mat> &input_imgs, const SfMData &input_sfmd, const Size &input_sz, const Mat** all_fundamental_matrices, vector<PolyLineGraph2DHMapImpl> &plgs, const vector<PolyLine2DMapSearch> &plmaps);
};

#endif /* INCLUDE_ALL_COMBINATIONS_CONSENSUS_MANAGER_PLGPCM3VIEWSPLGFOLLOWING_HPP_ */
