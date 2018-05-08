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

#include "segmented_edge_images_edge_manager.hpp"

#include "SfMData.h"
#include "convert_edge_images_pixel_to_segment.hpp"

using namespace std;
using namespace cv;
//: EdgeManager(vector<Mat> &imgs)
SegmentedEdgeImagesEdgeManager::SegmentedEdgeImagesEdgeManager(const vector<Mat> &imgs, const SfMData &input_sfmd, const cv::Mat** input_all_fundamental_matrices, const vector<Mat> &input_edge_images, const Vec<unsigned char, 3> &input_edge_color) : SegmentEdgeManager(imgs,input_sfmd,input_all_fundamental_matrices), edge_color(input_edge_color) {
		detect_edges(input_edge_images);
}

void SegmentedEdgeImagesEdgeManager::detect_edges(const vector<Mat> &imgs) {
	all_segments = convertEdgeImagePixelToSegments(imgs, edge_color);
}

SegmentedEdgeImagesEdgeManager::~SegmentedEdgeImagesEdgeManager() {}

