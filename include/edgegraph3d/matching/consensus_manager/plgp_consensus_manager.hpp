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


#ifndef INCLUDE_EDGEGRAPH3D_MATCHING_CONSENSUS_MANAGER_PLGP_CONSENSUS_MANAGER_HPP_
#define INCLUDE_EDGEGRAPH3D_MATCHING_CONSENSUS_MANAGER_PLGP_CONSENSUS_MANAGER_HPP_

#include <opencv2/core/mat.hpp>
#include <opencv2/core/types.hpp>
#include <utility>
#include <vector>

#include "polyline_graph_2d_hmap_impl.hpp"
#include "polyline_graph_2d.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"
#include "global_switches.hpp"
#include "glm.hpp"

using namespace cv;
using namespace std;

/**
 * Abstract parent class for all PLGPConsensusManager's implementations.
 *
 * A ConsensusManager addresses the task of determining a single 3D point from a set of
 * epipolar intersections that include outliers, by finding the 3D point most consistent
 * with the observations, if any. A PLGPConsensusManager is specialized for PLGPs.
 */
class PLGPConsensusManager {
public:
	/*
	 * return vector of new points tuples composed by
	 * 	(3d coordinates of the new point, vector of 2d reprojections on cameras where point is visible, ids of cameras where point is visible)
	 */
	virtual vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> consensus_strategy_single_point(const int starting_img_id, const int refpoint, const std::pair<std::vector<PolyLineGraph2D::plg_point>,std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>>> &intersections_and_correspondences_pair) =0;
	virtual vector<vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>>> consensus_strategy_single_point_vector(const int starting_img_id, const int refpoint, const std::pair<std::vector<PolyLineGraph2D::plg_point>,std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>>> &intersections_and_correspondences_pair) =0;
protected:
	const Mat** all_fundamental_matrices;
	vector<Mat> &imgs;
	const Size &img_size;
	const SfMData &sfmd;
	vector<PolyLineGraph2DHMapImpl> &plgs;
	PLGPConsensusManager(vector<Mat> &imgs, const SfMData &input_sfmd, const Size &input_sz, const Mat** all_fundamental_matrices, vector<PolyLineGraph2DHMapImpl> &plgs);
	virtual ~PLGPConsensusManager() {};
};

#endif /* INCLUDE_EDGE_MATCHER_CONSENSUS_MANAGER_CONSENSUSMANAGER_HPP_ */
