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


#ifndef INCLUDE_EDGEGRAPH3D_MATCHING_POLYLINE_MATCHING_POLYLINE_MATCHER_HPP_
#define INCLUDE_EDGEGRAPH3D_MATCHING_POLYLINE_MATCHING_POLYLINE_MATCHER_HPP_

#include <opencv2/core/types.hpp>
#include <set>
#include <utility>
#include <vector>

#include "polyline_graph_2d_hmap_impl.hpp"
#include "polyline_graph_2d.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"

#define FIND_WITHIN_DIST 10

pair<vector<ulong>,vector<vector<set<ulong>>>> polyline_matching_closeness_to_refpoints(const vector<PolyLineGraph2DHMapImpl> &plgs, const SfMData &sfmd, const Size &img_sz);

// Assumes close_refpoints_a and close_refpoints_b are both sorted
float compute_compatibility(const vector<float> &refpoints_weights, const vector<ulong> &close_refpoints_a_visible_on_img_b, const vector<ulong> &close_refpoints_b_visible_on_img_a);

float compute_refpoint_weight(const vector<set<ulong>> &close_polylines);

vector<long> read_cluster_info(const char *input_filename);

vector<vector<set<ulong>>> compute_polyline_matches_from_nodes_component_ids(const vector<pair<int,ulong>> &polyline_matches_vector,int amount_of_plgs,const vector<long> &nodes_component_ids);

/* Compute close components:
 * Output:
 *	- set of polyline_id for each close polyline to refpoint on img, for each img, for each refpoint
 *	- vector of close refpoints (sorted in ascending order), for each polyline, for each image
 *	- for each component in the potential compatibility graph (i.e. elements both close to at least a refpoint): for each image the potential matches (polyline ids)
 */
tuple<vector<vector<set<ulong>>>,vector<vector<vector<ulong>>>,vector<vector<set<ulong>>>> polyline_matching_similarity_graph(const vector<PolyLineGraph2DHMapImpl> &plgs, const SfMData &sfmd, const Size &img_sz, const char *out_compatibility_graph_file, const char *out_polyline_communities);

#endif /* INCLUDE_EDGEGRAPH3D_MATCHING_POLYLINE_MATCHING_POLYLINE_MATCHER_HPP_ */
