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

#include "edge_manager.hpp"
#include "glm.hpp"

std::vector<std::vector<std::vector<glm::vec2>>> EdgeManager::detect_epipolar_intersections_for_all_starting_intersections(const int starting_image_id, const int starting_point_id, const vector<glm::vec2> &intersection_points_on_starting_image, const float max_correspondence_detection_radius) {
	std::vector<std::vector<std::vector<glm::vec2>>> result;

	for (const auto &intersection_point_on_starting_image : intersection_points_on_starting_image)
		result.push_back(detect_epipolar_intersections(starting_image_id, starting_point_id, intersection_point_on_starting_image, max_correspondence_detection_radius));

	return result;
}

std::pair<std::vector<glm::vec2>,std::vector<std::vector<std::vector<glm::vec2>>>> EdgeManager::detect_nearby_intersections_and_correspondences(const int starting_image_id, const int starting_point_id, const float starting_detection_dist, const float max_correspondence_detection_radius) {
	std::vector<glm::vec2> nearby_intersections = detect_nearby_edge_intersections(starting_image_id, starting_point_id, starting_detection_dist);
	std::vector<std::vector<std::vector<glm::vec2>>> all_correspondences = detect_epipolar_intersections_for_all_starting_intersections(starting_image_id, starting_point_id, nearby_intersections, max_correspondence_detection_radius);
	return std::pair<std::vector<glm::vec2>,std::vector<std::vector<std::vector<glm::vec2>>>>(nearby_intersections,all_correspondences);
}


