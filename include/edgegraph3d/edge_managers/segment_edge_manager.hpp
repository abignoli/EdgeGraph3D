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


#ifndef INCLUDE_EDGEGRAPH3D_EDGE_MANAGERS_SEGMENT_EDGE_MANAGER_HPP_
#define INCLUDE_EDGEGRAPH3D_EDGE_MANAGERS_SEGMENT_EDGE_MANAGER_HPP_

#include <utility>
#include <vector>

#include "glm.hpp"

#include "opencv2/core/utility.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"

#include "edge_manager.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"

using namespace std;
using namespace cv;

// If enabled closest points (line,segment) while computing epipolar intersections if no actual intersection is found, up to MAX_CLOSE_POINT_DISTANCE
#define SEGMENT_EDGE_MANAGER_SELECT_CLOSE_POINTS_ENABLED

#define MAX_CLOSE_POINT_DISTANCE 1
#define MAX_CLOSE_POINT_DISTANCESQ (MAX_CLOSE_POINT_DISTANCE*MAX_CLOSE_POINT_DISTANCE)

class SegmentEdgeManager : public EdgeManager {
public:
	struct segment_point {
		glm::vec2 coords;
		glm::vec4 segment_extremes;

		segment_point() {}
		segment_point(glm::vec2 input_coords, glm::vec4 input_segment_extremes) : coords(input_coords), segment_extremes(input_segment_extremes) {}
	};

	std::vector<glm::vec2> detect_nearby_edge_intersections(const int imgId, const int startingPoint_id, const float starting_detection_dist);
	std::vector<std::vector<glm::vec2>> detect_epipolar_intersections(const int starting_image_id, const int starting_point_id, const glm::vec2 &intersection_point_on_starting_image, const float max_correspondence_detection_radius);
	std::vector<std::vector<std::vector<glm::vec2>>> detect_epipolar_intersections_for_all_starting_intersections(const int starting_image_id, const int starting_point_id, const vector<glm::vec2> &intersection_points_on_starting_image, const float max_correspondence_detection_radius);
	std::pair<std::vector<glm::vec2>,std::vector<std::vector<std::vector<glm::vec2>>>> detect_nearby_intersections_and_correspondences(const int starting_image_id, const int starting_point_id, const float starting_detection_dist, const float max_correspondence_detection_radius);

	std::vector<segment_point> detect_nearby_edge_intersections_with_segment_output(const int imgId, const int startingPoint_id, const float starting_detection_dist);
	std::vector<std::vector<segment_point>> detect_epipolar_intersections_with_segment_output(const int starting_image_id, const int starting_point_id, const segment_point &intersection_point_on_starting_image, const float max_correspondence_detection_radius);
	std::vector<std::vector<std::vector<segment_point>>> detect_epipolar_intersections_for_all_starting_intersections_with_segment_output(const int starting_image_id, const int starting_point_id, const vector<segment_point> &intersection_points_on_starting_image, const float max_correspondence_detection_radius);
	std::pair<std::vector<segment_point>,std::vector<std::vector<std::vector<segment_point>>>> detect_nearby_intersections_and_correspondences_with_segment_output(const int starting_image_id, const int starting_point_id, const float starting_detection_dist, const float max_correspondence_detection_radius);

	vector<Mat> get_edge_images();
	vector<Mat> get_edge_images_original_background();
	void find_closest_segment_projection(const int img_id, const glm::vec2 &point, glm::vec2 &closest_projection, float &min_distance, glm::vec4 &closest_segment);
protected:
	const vector<Mat> greyscale_imgs;
	vector<vector<glm::vec4>> all_segments; // all line segments for each image
	virtual void detect_edges(const vector<Mat> &imgs) =0;
	SegmentEdgeManager(const vector<Mat> &imgs, const SfMData &input_sfmd, const cv::Mat** input_all_fundamental_matrices);
	virtual ~SegmentEdgeManager() {};
private:
	std::vector<glm::vec2> detect_epipolar_intersections_on_image(const int img_to_find_intersections_on, const glm::vec2 &starting_point_on_current_image, const float max_correspondence_detection_radius, const glm::vec3 &epipolar_line);

	std::vector<segment_point> detect_epipolar_intersections_on_image_with_segment_output(const int img_to_find_intersections_on, const glm::vec2 &starting_point_on_current_image, const float max_correspondence_detection_radius, const glm::vec3 &epipolar_line);

};


#endif /* INCLUDE_EDGEGRAPH3D_EDGE_MANAGERS_SEGMENT_EDGE_MANAGER_HPP_ */
