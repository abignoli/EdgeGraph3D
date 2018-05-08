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

#include <opencv2/core/core.hpp>
#include <opencv2/core/cvdef.h>
#include <opencv2/core/cvstd.hpp>
#include <opencv2/core/mat.inl.hpp>
#include <opencv2/core/matx.hpp>
#include <opencv2/core/types.hpp>
#include <opencv2/imgproc/types_c.h>
#include <opencv2/imgproc.hpp>

#include <cmath>
#include <limits>

#include <glm.hpp>
#include "segment_edge_manager.hpp"

#include "SfMData.h"
#include "edge_graph_3d_utilities.hpp"
#include "geometric_utilities.hpp"

using namespace std;
using namespace cv;
//: EdgeManager(vector<Mat> &imgs)
SegmentEdgeManager::SegmentEdgeManager(const vector<Mat> &imgs, const SfMData &input_sfmd, const cv::Mat** input_all_fundamental_matrices) : EdgeManager(imgs,input_sfmd,input_all_fundamental_matrices)  {
	vector<Mat> gs;

	for (const Mat& img: imgs) {
		Mat cur_gs;
		cvtColor( img, cur_gs, CV_BGR2GRAY );
		gs.push_back(cur_gs);
	}

	greyscale_imgs = gs;
}

std::vector<std::vector<glm::vec2>> SegmentEdgeManager::detect_epipolar_intersections(const int starting_image_id, const int starting_point_id, const glm::vec2 &intersection_point_on_starting_image, const float max_correspondence_detection_radius) {
	std::vector<std::vector<glm::vec2>> all_epipolar_intersections;
	std::vector<glm::vec2> cur_epipolar_intersections;
	glm::vec2 starting_point_on_cur_img;
	glm::vec3 epipolar_line;
	int cams_observing_point_amount = sfmd.camViewingPointN_[starting_point_id].size();
	int cur_img;

	for(int i=0; i < cams_observing_point_amount; i++) {
		cur_img = sfmd.camViewingPointN_[starting_point_id][i];
		if (cur_img != starting_image_id) {
			// compute starting_point_on_cur_img
			starting_point_on_cur_img = sfmd.point2DoncamViewingPoint_[starting_point_id][i];

			// compute epipolar line
			const cv::Mat &F = all_fundamental_matrices[starting_image_id][cur_img];
			if(!computeCorrespondEpilineSinglePoint(intersection_point_on_starting_image, F, epipolar_line,1)) {
				all_epipolar_intersections.push_back(std::vector<glm::vec2>());
				continue;
			}

			// compute current image correspondences
			cur_epipolar_intersections = detect_epipolar_intersections_on_image(cur_img, starting_point_on_cur_img, max_correspondence_detection_radius, epipolar_line);

			// add current image correspondences
			all_epipolar_intersections.push_back(cur_epipolar_intersections);
		} else {
			// element corresponding to starting image in the defined order, just put a vector containing the starting intersection point coordinates
			std::vector<glm::vec2> original_intersection;
			original_intersection.push_back(intersection_point_on_starting_image);
			all_epipolar_intersections.push_back(original_intersection);
		}
	}

	return all_epipolar_intersections;
}

std::vector<glm::vec2> SegmentEdgeManager::detect_nearby_edge_intersections(const int imgId, const int startingPoint_id, const float starting_detection_dist) {
	std::vector<glm::vec2> nearby_intersections;
	std::vector<glm::vec2> cur_nearby_intersections;
	bool tmp;
	glm::vec2 startingPoint;
	const vector<glm::vec4> &cur_segments = all_segments[imgId];

	// get 2D coordinates of starting point on starting image
	get_2d_coordinates_of_point_on_image(sfmd, imgId, startingPoint_id, startingPoint, tmp);

	for(const auto &segm: cur_segments) {
		cur_nearby_intersections = detect_circle_segment_intersections(segm,startingPoint,starting_detection_dist);
		nearby_intersections.insert(std::end(nearby_intersections), std::begin(cur_nearby_intersections), std::end(cur_nearby_intersections));
	}

	return nearby_intersections;
}

std::vector<glm::vec2> SegmentEdgeManager::detect_epipolar_intersections_on_image(const int img_to_find_intersections_on, const glm::vec2 &starting_point_on_current_image, const float max_correspondence_detection_radius, const glm::vec3 &epipolar_line) {
	std::vector<glm::vec2> intersections;
	glm::vec2 cur_intersection;
	bool found;
	vector<glm::vec4> &cur_segments = all_segments[img_to_find_intersections_on];
	float detection_dist_sq = pow(max_correspondence_detection_radius,2);

	for(const auto &segm: cur_segments) {

#ifdef SEGMENT_EDGE_MANAGER_SELECT_CLOSE_POINTS_ENABLED
		single_intersect_segment_line_with_close_points_on_segm_in_detection_range(segm,epipolar_line, starting_point_on_current_image,detection_dist_sq, MAX_CLOSE_POINT_DISTANCESQ, cur_intersection, found);
#else
		single_intersect_segment_line_in_detection_range(segm,epipolar_line, starting_point_on_current_image,detection_dist_sq, cur_intersection,found);
#endif

		if (found)
			intersections.push_back(cur_intersection);
	}

	return intersections;
}

vector<Mat> SegmentEdgeManager::get_edge_images() {
	vector<Mat> edge_images;
	glm::vec4 e;

	Ptr<LineSegmentDetector> ls = createLineSegmentDetector(LSD_REFINE_ADV);

	for(int i=0; i< all_segments.size(); i++) {
		Mat edge_image = Mat(original_imgs[0].rows,original_imgs[0].cols,CV_8UC3);
		edge_image.setTo(cv::Scalar(0,0,0));
		vector<Vec4f> cur_edges = convert_vec4_Vec4f(all_segments[i]);
		//cout << "e " << cur_edges[0] << endl;
		ls->drawSegments(edge_image, cur_edges);
		edge_images.push_back(edge_image);
	}
	return edge_images;
}

vector<Mat> SegmentEdgeManager::get_edge_images_original_background() {
	vector<Mat> edge_images;
	glm::vec4 e;

	Ptr<LineSegmentDetector> ls = createLineSegmentDetector(LSD_REFINE_ADV);

	for(int i=0; i< all_segments.size(); i++) {
		Mat edge_image = greyscale_imgs[i].clone();
		vector<Vec4f> cur_edges = convert_vec4_Vec4f(all_segments[i]);
		//cout << "e " << cur_edges[0] << endl;
		ls->drawSegments(edge_image, cur_edges);
		edge_images.push_back(edge_image);
	}

	return edge_images;
}

/**
 * Given an image and a 2D point on it, find the closest segment along with the closest projection of the point on the segment and its distance.
 */
void SegmentEdgeManager::find_closest_segment_projection(const int img_id, const glm::vec2 &point, glm::vec2 &closest_projection, float &min_distance, glm::vec4 &closest_segment) {
	vector<glm::vec4> &img_segments = all_segments[img_id];
	glm::vec2 cur_projection;
	float dist;

	min_distance=std::numeric_limits<float>::max();

	for(const auto &segm : img_segments) {
		dist = minimum_distancesq(segm,point, cur_projection);
		if(dist<min_distance) {
			closest_projection = cur_projection;
			closest_segment = segm;
			min_distance = dist;
		}
	}

	min_distance = sqrt(min_distance);
}






std::vector<SegmentEdgeManager::segment_point> SegmentEdgeManager::detect_epipolar_intersections_on_image_with_segment_output(const int img_to_find_intersections_on, const glm::vec2 &starting_point_on_current_image, const float max_correspondence_detection_radius, const glm::vec3 &epipolar_line) {
	std::vector<SegmentEdgeManager::segment_point> intersections;
	glm::vec2 cur_intersection;
	bool found;
	vector<glm::vec4> &cur_segments = all_segments[img_to_find_intersections_on];
	float detection_dist_sq = pow(max_correspondence_detection_radius,2);

	for(const auto &segm: cur_segments) {

#ifdef SELECT_CLOSE_POINTS_ENABLED
		single_intersect_segment_line_with_close_points_on_segm_in_detection_range(segm,epipolar_line, starting_point_on_current_image,detection_dist_sq, MAX_CLOSE_POINT_DISTANCESQ, cur_intersection, found);
#else
		single_intersect_segment_line_in_detection_range(segm,epipolar_line, starting_point_on_current_image,detection_dist_sq, cur_intersection,found);
#endif

		if (found)
			intersections.push_back(SegmentEdgeManager::segment_point(cur_intersection,segm));

	}

	return intersections;
}

std::vector<SegmentEdgeManager::segment_point> SegmentEdgeManager::detect_nearby_edge_intersections_with_segment_output(const int imgId, const int startingPoint_id, const float starting_detection_dist) {
	std::vector<SegmentEdgeManager::segment_point> nearby_intersections;
	std::vector<glm::vec2> cur_nearby_intersections;
	bool tmp;
	glm::vec2 startingPoint;
	const vector<glm::vec4> &cur_segments = all_segments[imgId];

	// get 2D coordinates of starting point on starting image
	get_2d_coordinates_of_point_on_image(sfmd, imgId, startingPoint_id, startingPoint, tmp);

	for(const auto &segm: cur_segments) {
		cur_nearby_intersections = detect_circle_segment_intersections(segm,startingPoint,starting_detection_dist);

		for(const auto &intersection : cur_nearby_intersections)
			nearby_intersections.push_back(SegmentEdgeManager::segment_point(intersection,segm));
	}

	return nearby_intersections;
}

std::vector<std::vector<SegmentEdgeManager::segment_point>> SegmentEdgeManager::detect_epipolar_intersections_with_segment_output(const int starting_image_id, const int starting_point_id, const SegmentEdgeManager::segment_point &intersection_point_on_starting_image, const float max_correspondence_detection_radius) {
	glm::vec2 intersection_point_on_starting_image_coords = intersection_point_on_starting_image.coords;
	std::vector<std::vector<SegmentEdgeManager::segment_point>> all_epipolar_intersections;
	std::vector<SegmentEdgeManager::segment_point> cur_epipolar_intersections;
	glm::vec2 starting_point_on_cur_img;
	glm::vec3 epipolar_line;
	int cams_observing_point_amount = sfmd.camViewingPointN_[starting_point_id].size();
	int cur_img;

	for(int i=0; i < cams_observing_point_amount; i++) {
		cur_img = sfmd.camViewingPointN_[starting_point_id][i];
		if (cur_img != starting_image_id) {
			// compute starting_point_on_cur_img
			starting_point_on_cur_img = sfmd.point2DoncamViewingPoint_[starting_point_id][i];

			// compute epipolar line
			const cv::Mat &F = all_fundamental_matrices[starting_image_id][cur_img];

			if(!computeCorrespondEpilineSinglePoint(intersection_point_on_starting_image_coords, F, epipolar_line,1)) {
				all_epipolar_intersections.push_back(std::vector<SegmentEdgeManager::segment_point>());
				continue;
			}
			// compute current image correspondences
			cur_epipolar_intersections = detect_epipolar_intersections_on_image_with_segment_output(cur_img, starting_point_on_cur_img, max_correspondence_detection_radius, epipolar_line);

			// add current image correspondences
			all_epipolar_intersections.push_back(cur_epipolar_intersections);
		} else {
			// element corresponding to starting image in the defined order, just put a vector containing the starting intersection point coordinates
			std::vector<SegmentEdgeManager::segment_point> original_intersection;
			original_intersection.push_back(intersection_point_on_starting_image);
			all_epipolar_intersections.push_back(original_intersection);
		}
	}

	return all_epipolar_intersections;
}

std::vector<std::vector<std::vector<SegmentEdgeManager::segment_point>>> SegmentEdgeManager::detect_epipolar_intersections_for_all_starting_intersections_with_segment_output(const int starting_image_id, const int starting_point_id, const vector<SegmentEdgeManager::segment_point> &intersection_points_on_starting_image, const float max_correspondence_detection_radius) {
	std::vector<std::vector<std::vector<SegmentEdgeManager::segment_point>>> result;

	for (const auto &intersection_point_on_starting_image : intersection_points_on_starting_image)
		result.push_back(detect_epipolar_intersections_with_segment_output(starting_image_id, starting_point_id, intersection_point_on_starting_image, max_correspondence_detection_radius));

	return result;
}

std::pair<std::vector<SegmentEdgeManager::segment_point>,std::vector<std::vector<std::vector<SegmentEdgeManager::segment_point>>>> SegmentEdgeManager::detect_nearby_intersections_and_correspondences_with_segment_output(const int starting_image_id, const int starting_point_id, const float starting_detection_dist, const float max_correspondence_detection_radius) {
	std::vector<SegmentEdgeManager::segment_point> nearby_intersections = detect_nearby_edge_intersections_with_segment_output(starting_image_id, starting_point_id, starting_detection_dist);
	std::vector<std::vector<std::vector<SegmentEdgeManager::segment_point>>> all_correspondences = detect_epipolar_intersections_for_all_starting_intersections_with_segment_output(starting_image_id, starting_point_id, nearby_intersections, max_correspondence_detection_radius);
	return std::pair<std::vector<SegmentEdgeManager::segment_point>,std::vector<std::vector<std::vector<SegmentEdgeManager::segment_point>>>>(nearby_intersections,all_correspondences);
}
