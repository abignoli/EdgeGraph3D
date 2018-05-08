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


#ifndef EDGEMANAGER_H_
#define EDGEMANAGER_H_

#include <opencv2/core/mat.hpp>
#include <utility>
#include <vector>

#include "global_defines.hpp"
#include "global_switches.hpp"
#include "glm.hpp"


struct SfMData;

using namespace std;
using namespace cv;

/**
 * Abstract parent class for all EdgeManager's implementations.
 *
 * An EdgeManager has the task of managing operations related to edges on the input dataset,
 * such as edge detection and edge intersections detection.
 */
class EdgeManager {
public:
	virtual std::vector<glm::vec2> detect_nearby_edge_intersections(const int imgId, const int startingPoint_id, const float starting_detection_dist) =0;
	virtual std::vector<std::vector<glm::vec2>> detect_epipolar_intersections(const int starting_image_id, const int starting_point_id, const glm::vec2 &intersection_point_on_starting_image, const float max_correspondence_detection_radius) =0;
	std::vector<std::vector<std::vector<glm::vec2>>> detect_epipolar_intersections_for_all_starting_intersections(const int starting_image_id, const int starting_point_id, const vector<glm::vec2> &intersection_points_on_starting_image, const float max_correspondence_detection_radius);
	std::pair<std::vector<glm::vec2>,std::vector<std::vector<std::vector<glm::vec2>>>> detect_nearby_intersections_and_correspondences(const int starting_image_id, const int starting_point_id, const float starting_detection_dist, const float max_correspondence_detection_radius);
	virtual vector<Mat> get_edge_images() =0;
	virtual vector<Mat> get_edge_images_original_background() =0;
protected:
	const vector<Mat> original_imgs;
	const SfMData &sfmd;
	const cv::Mat** all_fundamental_matrices;

	virtual void detect_edges(const vector<Mat> &imgs) =0;
	virtual std::vector<glm::vec2> detect_epipolar_intersections_on_image(const int img_to_find_intersections_on, const glm::vec2 &starting_point_on_current_image, const float max_correspondence_detection_radius, const glm::vec3 &epipolar_line) =0;
	EdgeManager(const vector<Mat> &imgs, const SfMData &input_sfmd, const cv::Mat** input_all_fundamental_matrices) : original_imgs(imgs), sfmd(input_sfmd),all_fundamental_matrices(input_all_fundamental_matrices) {
		//cout << "EdgeManager constructor call\n";
	}
	virtual ~EdgeManager() {};
};

#endif /* EDGEMANAGER_H_ */
