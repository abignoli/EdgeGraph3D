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


#ifndef INCLUDE_EDGEGRAPH3D_EDGE_MANAGERS_PLG_EDGE_MANAGER_HPP_
#define INCLUDE_EDGEGRAPH3D_EDGE_MANAGERS_PLG_EDGE_MANAGER_HPP_

#include <utility>
#include <vector>

#include "glm.hpp"

#include <opencv2/core/mat.hpp>

#include "polyline_graph_2d_hmap_impl.hpp"
#include "opencv2/core/utility.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"

#include "edge_manager.hpp"
#include "global_defines.hpp"
#include "polyline_2d_map_search.hpp"
#include "polyline_graph_2d.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"
#include "global_switches.hpp"

using namespace std;
using namespace cv;

class PLGEdgeManager : public EdgeManager {
public:

	std::vector<glm::vec2> detect_nearby_edge_intersections(const int imgId, const int startingPoint_id, const float starting_detection_dist);
	std::vector<std::vector<glm::vec2>> detect_epipolar_intersections(const int starting_image_id, const int starting_point_id, const glm::vec2 &intersection_point_on_starting_image, const float max_correspondence_detection_radius);
	std::vector<std::vector<std::vector<glm::vec2>>> detect_epipolar_intersections_for_all_starting_intersections(const int starting_image_id, const int starting_point_id, const vector<glm::vec2> &intersection_points_on_starting_image, const float max_correspondence_detection_radius);
	std::pair<std::vector<glm::vec2>,std::vector<std::vector<std::vector<glm::vec2>>>> detect_nearby_intersections_and_correspondences(const int starting_image_id, const int starting_point_id, const float starting_detection_dist, const float max_correspondence_detection_radius);


	std::vector<PolyLineGraph2D::plg_point> detect_nearby_edge_intersections_plgp(const int imgId, const int startingPoint_id);
	std::vector<std::vector<PolyLineGraph2D::plg_point>> detect_epipolar_intersections_plgp(const int starting_image_id, const int starting_point_id, const PolyLineGraph2D::plg_point &intersection_point_on_starting_image_plgp, const float max_correspondence_detection_radius);
	std::vector<PolyLineGraph2D::plg_point> detect_epipolar_intersections_on_image_plgp(const int img_to_find_intersections_on, const glm::vec2 &starting_point_on_current_image, const float max_correspondence_detection_radius, const glm::vec3 &epipolar_line);
	std::pair<std::vector<PolyLineGraph2D::plg_point>,std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>>> detect_nearby_intersections_and_correspondences_plgp(const int starting_image_id, const int starting_point_id);

	std::vector<PolyLineGraph2D::plg_point> detect_epipolar_intersections_on_image_plgp(const int img_to_find_intersections_on, const glm::vec2 &starting_point_on_current_image, const float max_correspondence_detection_radius, const glm::vec3 &epipolar_line, const vector<ulong> &potentially_correspondent_polylines);
	std::vector<std::vector<PolyLineGraph2D::plg_point>> detect_epipolar_intersections_plgp(const int starting_image_id, const int starting_point_id, const PolyLineGraph2D::plg_point &intersection_point_on_starting_image_plgp, const float max_correspondence_detection_radius, const vector<vector<ulong>> &potentially_correspondent_polylines);
	std::vector<std::pair<std::vector<PolyLineGraph2D::plg_point>,std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>>>> detect_nearby_intersections_and_correspondences_plgp(const int starting_point_id);

	std::vector<std::pair<std::vector<PolyLineGraph2D::plg_point>,std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>>>> detect_nearby_intersections_and_correspondences_plgp_closest_starting_plgp(const int starting_point_id);

	vector<Mat> get_edge_images();
	vector<Mat> get_edge_images_original_background();
	PLGEdgeManager(const vector<Mat> &imgs, const SfMData &input_sfmd, const cv::Mat** input_all_fundamental_matrices, vector<PolyLineGraph2DHMapImpl> &input_plgs, const float starting_detection_dist, const float correspondence_detection_range_multiplication_factor);
	~PLGEdgeManager() {};
protected:
	//vector<PolyLine2DMapSearch> starting_plmaps;
	vector<PolyLine2DMapSearch> correspondence_plmaps;
	const float starting_detection_dist;
	const float starting_detection_distsq;
	const float correspondence_detection_range_multiplication_factor;
	const float correspondence_detection_dist;
	const float correspondence_detection_distsq;
	vector<PolyLineGraph2DHMapImpl> &plgs;
	const vector<Mat> greyscale_imgs;
	vector<vector<glm::vec4>> all_segments; // all line segments for each image
	vector<vector<pair<ulong,ulong>>> all_segments_polyline_ids; // for each image, for each segment, its (polyline_id, segment_index)
	void detect_edges(const vector<Mat> &imgs);
private:
	std::vector<glm::vec2> detect_epipolar_intersections_on_image(const int img_to_find_intersections_on, const glm::vec2 &starting_point_on_current_image, const float max_correspondence_detection_radius, const glm::vec3 &epipolar_line);

	std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>> detect_epipolar_intersections_for_all_starting_intersections_plgp(const int starting_image_id, const int starting_point_id, const vector<PolyLineGraph2D::plg_point> &intersection_points_on_starting_image);

	std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>> detect_epipolar_intersections_for_all_starting_intersections_plgp(const int starting_image_id, const int starting_point_id, const vector<PolyLineGraph2D::plg_point> &intersection_points_on_starting_image, const vector<vector<ulong>> &potentially_correspondent_polylines);
};



#endif /* INCLUDE_EDGE_MATCHER_EDGE_MANAGER_PLGEDGEMANAGERCLOSESTONLYCLOSESTONLY_HPP_ */
