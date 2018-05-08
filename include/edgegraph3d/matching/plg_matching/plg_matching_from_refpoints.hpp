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


#ifndef INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_PLG_MATCHING_FROM_REFPOINTS_HPP_
#define INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_PLG_MATCHING_FROM_REFPOINTS_HPP_

#include <opencv2/core/mat.hpp>
#include <vector>

#include "polyline_graph_2d.hpp"
#include "polyline_graph_2d_hmap_impl.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"
#include "global_switches.hpp"
#include "glm.hpp"

class PLGEdgeManagerClosestOnly;
class PLGMatchesManager;
class PLGPConsensusManager;

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> plg_matching_from_refpoint_starting_image(const SfMData &sfm_data, const EdgeManager *em,const PLGPConsensusManager *cm, const ulong refpoint_id, const ulong starting_img_id);

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> plg_matching_from_refpoint(const SfMData &sfm_data, const EdgeManager *em,const PLGPConsensusManager *cm, const ulong refpoint_id, PLGMatchesManager &plgmm);

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> plg_matching_from_refpoints_parallel(const SfMData &sfm_data, const EdgeManager *em,const PLGPConsensusManager *cm, PLGMatchesManager &plgmm);

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> plg_matching_from_refpoints(const SfMData &sfm_data, const EdgeManager *em,const PLGPConsensusManager *cm, PLGMatchesManager &plgmm);

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> plg_matching_from_refpoint(const SfMData &sfm_data, const EdgeManager *em,const PLGPConsensusManager *cm, const ulong refpoint_id);

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> plg_matching_from_refpoints_parallel(const SfMData &sfm_data, const EdgeManager *em,const PLGPConsensusManager *cm);

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> plg_matching_from_refpoints(const SfMData &sfm_data, const EdgeManager *em,const PLGPConsensusManager *cm);

vector<std::tuple<glm::vec3,vector<glm::vec2>,vector<int>>> process_refpoint(const vector<Mat> &imgs, vector<Mat> &out_imgs, const vector<PolyLineGraph2DHMapImpl> &plgs, const SfMData &sfmd, const Mat** all_fundamental_matrices, const PLGEdgeManagerClosestOnly *em, const ulong refpoint_id, char * outfolder);

#endif /* INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_PLG_MATCHING_FROM_REFPOINTS_HPP_ */
