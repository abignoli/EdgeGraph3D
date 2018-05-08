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


#include "plg_analysis_utilities.hpp"

#include <opencv2/core/cvdef.h>
#include <opencv2/core/mat.hpp>
#include <opencv2/core/mat.inl.hpp>
#include <opencv2/core/types.hpp>
#include <opencv2/imgcodecs.hpp>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>

#include "SfMData.h"
#include "convert_edge_images_pixel_to_segment.hpp"
#include "glm.hpp"

#include "polyline_graph_2d.hpp"
#include "global_defines.hpp"
#include "drawing_utilities.hpp"
#include "geometric_utilities.hpp"
#include "edge_graph_3d_utilities.hpp"

int is_edgerefpoint(const vector<vector<glm::vec4>> &vec_segments, const SfMData &sfmd, int refpoint_id) {
	glm::vec2 projection;

	for(int i=0; i < sfmd.camViewingPointN_[refpoint_id].size(); i++) {
		int img_id = sfmd.camViewingPointN_[refpoint_id][i];
		const glm::vec2 &p_coords = sfmd.point2DoncamViewingPoint_[refpoint_id][i];
		vector<glm::vec4> segments = vec_segments[img_id];
		for(const auto &s : segments)
			if(minimum_distancesq(s, p_coords, projection) < MAX_EDGEREFPOINT_DISTSQ)
				return true;
	}

	return false;
}

std::vector<int> find_edgerefpoints(const vector<PolyLineGraph2DHMapImpl> &plgs, const SfMData &sfmd) {
	std::vector<int> res;

	vector<vector<glm::vec4>> vec_segments;
	for(const auto &plg : plgs)
		vec_segments.push_back(plg.get_segments_list());

	for(int i=0; i < sfmd.numPoints_; i++)
		if(is_edgerefpoint(vec_segments, sfmd, i))
			res.push_back(i);

	return res;
}
