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


#ifndef INCLUDE_EDGEGRAPH3D_EDGE_MANAGERS_INPUT_SEGMENTS_EDGE_MANAGER_HPP_
#define INCLUDE_EDGEGRAPH3D_EDGE_MANAGERS_INPUT_SEGMENTS_EDGE_MANAGER_HPP_


#include <opencv2/core/mat.hpp>
#include <opencv2/core/matx.hpp>
#include <vector>

#include "glm.hpp"
#include "opencv2/core/utility.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"

#include "segment_edge_manager.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"

using namespace std;
using namespace cv;

class InputSegmentsEdgeManager : public SegmentEdgeManager {
public:
	InputSegmentsEdgeManager(const vector<Mat> &imgs, const SfMData &input_sfmd, const cv::Mat** input_all_fundamental_matrices, vector<vector<glm::vec4>> input_vec_segments);
	~InputSegmentsEdgeManager();
private:
	Vec<unsigned char, 3> edge_color;
	void detect_edges(const vector<Mat> &imgs);
};


#endif /* INCLUDE_EDGEGRAPH3D_EDGE_MANAGERS_INPUT_SEGMENTS_EDGE_MANAGER_HPP_ */
