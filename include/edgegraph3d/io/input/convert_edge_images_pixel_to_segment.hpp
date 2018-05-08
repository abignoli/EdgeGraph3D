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


#ifndef INCLUDE_EDGEGRAPH3D_IO_INPUT_CONVERT_EDGE_IMAGES_PIXEL_TO_SEGMENT_HPP_
#define INCLUDE_EDGEGRAPH3D_IO_INPUT_CONVERT_EDGE_IMAGES_PIXEL_TO_SEGMENT_HPP_

#include <opencv2/core/mat.hpp>
#include <opencv2/core/matx.hpp>
#include <vector>

#include "polyline_graph_2d_hmap_impl.hpp"
#include "graph_adjacency_list_undirected.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"
#include "global_switches.hpp"
#include "glm.hpp"

using namespace std;
using namespace cv;

/**
 * Takes an edge image, assigns a Node ID to each edge pixel in the image
 */
ulong** convertEdgeImagesPixelToNodes(const Mat &img, const Vec<unsigned char, 3> &edge_color, ulong &amount_of_nodes);

GraphAdjacencyListUndirected<glm::vec2> convertEdgeImagePixelToGraph(const Mat &img, const Vec<unsigned char, 3> &edge_color);

vector<glm::vec4> convert_graph_to_edge_list(const GraphAdjacencyListUndirected<glm::vec2> &g);

vector<glm::vec4> convertEdgeImagePixelToSegments(const Mat &img, const Vec<unsigned char, 3> &edge_color);

Mat convertEdgeImagePixelToSegmentsImage(const Mat &img, const Vec<unsigned char, 3> &edge_color);

vector<Mat> convertEdgeImagesPixelToSegmentsImages(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color);

vector<GraphAdjacencyListUndirected<glm::vec2>> convertEdgeImagesPixelToGraphs(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color);

vector<vector<glm::vec4>> convertEdgeImagePixelToSegments(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color);

vector<Mat> convertEdgeImagesPixelToSegmentsImages_NoCycles(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color);

vector<Mat> convertEdgeImagesPixelToSegmentsImages_PolyLineGraph_simplified(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color);

vector<Mat> convertEdgeImagesPixelToSegmentsImages_MultiColor_PolyLineGraph_simplified(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color);

vector<Mat> convertEdgeImagesPixelToSegmentsImages_PolyLineGraph(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color);

vector<Mat> convertEdgeImagesPixelToSegmentsImages_MultiColor_PolyLineGraph(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color);

vector<Mat> convertEdgeImagesPixelToSegmentsImages_MultiColorPolyLines_PolyLineGraph_simplified(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color);

vector<Mat> convertEdgeImagesPixelToSegmentsImages_MultiColorComponents_PolyLineGraph_simplified(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color);

PolyLineGraph2DHMapImpl convertEdgeImagePolyLineGraph_simplified(const Mat &img, const Vec<unsigned char, 3> &edge_color);

vector<PolyLineGraph2DHMapImpl> convertEdgeImagesPolyLineGraphs_simplified(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color);

Mat convertEdgeImagePixelToSegmentsImage_MultiColorPolyLines_PolyLineGraph_simplified(const Mat &img, const Vec<unsigned char, 3> &edge_color);

vector<vector<glm::vec4>> convertEdgeImagePixelToSegments_PolyLineGraph_simplified(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color);

Mat convertEdgeImagePixelToSegmentsImage_MultiColorComponents_PolyLineGraph_optimized(const Mat &img, const Vec<unsigned char, 3> &edge_color);

vector<Mat> convertEdgeImagesPixelToSegmentsImages_MultiColorComponents_PolyLineGraph_optimized(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color);

/**
 * Adds:
 * - polyline simplification
 * - direct extreme connect
 */
PolyLineGraph2DHMapImpl convertEdgeImagePolyLineGraph_optimized(const Mat &img, const Vec<unsigned char, 3> &edge_color);

vector<PolyLineGraph2DHMapImpl> convert_edge_images_to_optimized_polyline_graphs(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color);

#endif /* INCLUDE_EDGEGRAPH3D_IO_INPUT_CONVERT_EDGE_IMAGES_PIXEL_TO_SEGMENT_HPP_ */
