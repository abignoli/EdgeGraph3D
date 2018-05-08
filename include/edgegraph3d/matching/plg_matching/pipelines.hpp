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


#ifndef INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_PIPELINES_HPP_
#define INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_PIPELINES_HPP_

#include <opencv2/core/mat.hpp>
#include <vector>

#include "polyline_graph_2d_hmap_impl.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"
#include "polyline_2d_map_search.hpp"

struct data_bundle;
struct edge_matcher_input_params;
class PLGMatchesManager;

using namespace std;
using namespace cv;

void edge_reconstruction_pipeline(const edge_matcher_input_params &emip, vector<Mat> &imgs, const char *out_folder, const vector<PolyLineGraph2DHMapImpl> &plgs, SfMData &sfmd, data_bundle *mfc, PLGMatchesManager &plgmm, const vector<PolyLine2DMapSearch> &plmaps);

#endif /* INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_PIPELINES_HPP_ */
