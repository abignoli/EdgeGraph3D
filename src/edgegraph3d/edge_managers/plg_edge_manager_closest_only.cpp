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
#include <opencv2/core/mat.inl.hpp>
#include <opencv2/imgproc/types_c.h>
#include <opencv2/imgproc.hpp>
#include <algorithm>
#include <cmath>

#include <glm.hpp>
#include "plg_edge_manager_closest_only.hpp"

#include "SfMData.h"
#include "drawing_utilities.hpp"
#include "geometric_utilities.hpp"
#include "NotImplementedException.hpp"

PLGEdgeManagerClosestOnly::PLGEdgeManagerClosestOnly(const vector<Mat> &imgs, const SfMData &input_sfmd, const cv::Mat** input_all_fundamental_matrices, vector<PolyLineGraph2DHMapImpl> &input_plgs) : EdgeManager(imgs,input_sfmd,input_all_fundamental_matrices), plgs(input_plgs)  {
	vector<Mat> gs;

	for (const Mat& img: imgs) {
		Mat cur_gs;
		cvtColor( img, cur_gs, CV_BGR2GRAY );
		gs.push_back(cur_gs);
	}

	greyscale_imgs = gs;

	for(const auto &plg : plgs) {
		vector<glm::vec4> segments_on_img;
		vector<pair<ulong,ulong>> polyline_ids_on_img;
		vector<pair<ulong,vector<glm::vec4>>> cur_img_segments = plg.get_segments_grouped_by_polyline_with_polyline_ids();
		for(const auto &polyline_pkg : cur_img_segments) {
			const ulong polyline_id = polyline_pkg.first;
			const vector<glm::vec4> &segments_on_polyline = polyline_pkg.second;
			for(ulong i=0; i < segments_on_polyline.size(); i++) {
				segments_on_img.push_back(segments_on_polyline[i]);
				polyline_ids_on_img.push_back(make_pair(polyline_id,i));
			}
		}
		all_segments.push_back(segments_on_img);
		all_segments_polyline_ids.push_back(polyline_ids_on_img);
	}
}

void PLGEdgeManagerClosestOnly::detect_edges(const vector<Mat> &imgs) {}

std::vector<std::vector<PolyLineGraph2D::plg_point>> PLGEdgeManagerClosestOnly::detect_epipolar_intersections_plgp(const int starting_image_id, const int starting_point_id, const PolyLineGraph2D::plg_point &intersection_point_on_starting_image_plgp, const float max_correspondence_detection_radius) {
	const glm::vec2 &intersection_point_on_starting_image = intersection_point_on_starting_image_plgp.plp.coords;
	std::vector<std::vector<PolyLineGraph2D::plg_point>> all_epipolar_intersections;
	std::vector<PolyLineGraph2D::plg_point> cur_epipolar_intersections;
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
				all_epipolar_intersections.push_back(std::vector<PolyLineGraph2D::plg_point>());
				continue;
			}
			// compute current image correspondences
			cur_epipolar_intersections = detect_epipolar_intersections_on_image_plgp(cur_img, starting_point_on_cur_img, max_correspondence_detection_radius, epipolar_line);

			// add current image correspondences
			all_epipolar_intersections.push_back(cur_epipolar_intersections);
		} else {
			// element corresponding to starting image in the defined order, just put a vector containing the starting intersection point coordinates
			std::vector<PolyLineGraph2D::plg_point> original_intersection;
			original_intersection.push_back(intersection_point_on_starting_image_plgp);
			all_epipolar_intersections.push_back(original_intersection);
		}
	}

	return all_epipolar_intersections;
}

std::vector<PolyLineGraph2D::plg_point> PLGEdgeManagerClosestOnly::detect_nearby_edge_intersections_plgp(const int imgId, const int startingPoint_id, const float starting_detection_dist) {
	// get only closest polyline if close enough

	glm::vec2 startingPoint;
	bool tmp;

	// get 2D coordinates of starting point on starting image
	get_2d_coordinates_of_point_on_image(sfmd, imgId, startingPoint_id, startingPoint, tmp);

	std::vector<PolyLineGraph2D::plg_point> res;
	PolyLineGraph2D::plg_point plgp;
	float min_dist;
	min_dist = plgs[imgId].cpf_find_unbound(startingPoint,plgp);

	if(min_dist < starting_detection_dist)
		 res.push_back(plgp);

	return res;
}

std::vector<PolyLineGraph2D::plg_point> PLGEdgeManagerClosestOnly::detect_epipolar_intersections_on_image_plgp(const int img_to_find_intersections_on, const glm::vec2 &starting_point_on_current_image, const float max_correspondence_detection_radius, const glm::vec3 &epipolar_line) {
	std::vector<PolyLineGraph2D::plg_point> intersections;
	glm::vec2 cur_intersection;
	bool found;
	const vector<glm::vec4> &cur_segments = all_segments[img_to_find_intersections_on];
	float detection_dist_sq = pow(max_correspondence_detection_radius,2);

	for(ulong i=0; i < cur_segments.size(); i++) {
		const glm::vec4 &segm = cur_segments[i];

#ifdef SELECT_CLOSE_POINTS_ENABLED
		single_intersect_segment_line_with_close_points_on_segm_in_detection_range(segm,epipolar_line, starting_point_on_current_image,detection_dist_sq, MAX_CLOSE_POINT_DISTANCESQ, cur_intersection, found);
#else
		single_intersect_segment_line_in_detection_range(segm,epipolar_line, starting_point_on_current_image,detection_dist_sq, cur_intersection,found);
#endif

		if (found)
			intersections.push_back(PolyLineGraph2D::plg_point(all_segments_polyline_ids[img_to_find_intersections_on][i],cur_intersection));

	}

	return intersections;
}

vector<Mat> PLGEdgeManagerClosestOnly::get_edge_images() {
	vector<Mat> edge_images;

	for(int i=0; i< all_segments.size(); i++) {
		Mat edge_image = draw_MultiColorComponents_PolyLineGraph_simplified(original_imgs[i], plgs[i]);
		edge_images.push_back(edge_image);
	}
	return edge_images;
}

vector<Mat> PLGEdgeManagerClosestOnly::get_edge_images_original_background() {
	vector<Mat> edge_images;

	for(int i=0; i< all_segments.size(); i++) {
		Mat edge_image = Mat(original_imgs[0].rows,original_imgs[0].cols,CV_8UC3);
		draw_overlay_MultiColorComponents_PolyLineGraph_simplified(edge_image, plgs[i]);
		edge_images.push_back(edge_image);
	}
	return edge_images;
}

std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>> PLGEdgeManagerClosestOnly::detect_epipolar_intersections_for_all_starting_intersections_plgp(const int starting_image_id, const int starting_point_id, const vector<PolyLineGraph2D::plg_point> &intersection_points_on_starting_image, const float correspondence_detection_range_multiplication_factor) {
	std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>> result;

	bool tmp;
	glm::vec2 init_coords;
	get_2d_coordinates_of_point_on_image(sfmd,starting_image_id,starting_point_id,init_coords,tmp);

	for (const auto &intersection_point_on_starting_image : intersection_points_on_starting_image) {
		float max_correspondence_detection_radius = compute_2d_distance(init_coords,intersection_point_on_starting_image.plp.coords) * correspondence_detection_range_multiplication_factor;
		result.push_back(detect_epipolar_intersections_plgp(starting_image_id, starting_point_id, intersection_point_on_starting_image, max_correspondence_detection_radius));
	}

	return result;
}

std::pair<std::vector<PolyLineGraph2D::plg_point>,std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>>> PLGEdgeManagerClosestOnly::detect_nearby_intersections_and_correspondences_plgp(const int starting_image_id, const int starting_point_id, const float starting_detection_dist, const float correspondence_detection_range_multiplication_factor) {
	std::vector<PolyLineGraph2D::plg_point> nearby_intersections = detect_nearby_edge_intersections_plgp(starting_image_id, starting_point_id, starting_detection_dist);
	std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>> all_correspondences = detect_epipolar_intersections_for_all_starting_intersections_plgp(starting_image_id, starting_point_id, nearby_intersections, correspondence_detection_range_multiplication_factor);
	return make_pair(nearby_intersections,all_correspondences);
}

void PLGEdgeManagerClosestOnly::detect_epipolar_intersections_plgp_exclude_parallel_epipsegments(const int starting_image_id, const int starting_point_id, const PolyLineGraph2D::plg_point &intersection_point_on_starting_image_plgp, const float max_correspondence_detection_radius, const float max_allowed_epip_abscos, std::vector<std::vector<PolyLineGraph2D::plg_point>> &res, bool &valid) {
	const glm::vec2 &intersection_point_on_starting_image = intersection_point_on_starting_image_plgp.plp.coords;
	std::vector<PolyLineGraph2D::plg_point> cur_epipolar_intersections;
	glm::vec2 starting_point_on_cur_img;
	glm::vec3 epipolar_line;
	int cams_observing_point_amount = sfmd.camViewingPointN_[starting_point_id].size();
	int cur_img;
	float cur_min_abscos;
	valid = true;

	res.clear();

	for(int i=0; i < cams_observing_point_amount; i++) {
		cur_img = sfmd.camViewingPointN_[starting_point_id][i];
		if (cur_img != starting_image_id) {
			// compute starting_point_on_cur_img
			starting_point_on_cur_img = sfmd.point2DoncamViewingPoint_[starting_point_id][i];

			// compute epipolar line
			const cv::Mat &F = all_fundamental_matrices[starting_image_id][cur_img];
			if(!computeCorrespondEpilineSinglePoint(intersection_point_on_starting_image, F, epipolar_line,1)) {
				res.push_back(std::vector<PolyLineGraph2D::plg_point>());
				continue;
			}
			// compute current image correspondences
			cur_epipolar_intersections = detect_epipolar_intersections_on_image_plgp_exclude_parallel_epipsegments(cur_img, starting_point_on_cur_img, max_correspondence_detection_radius, epipolar_line, cur_min_abscos);

			// If an intersection with almost parallel epipolar is found, interrupt
			if(cur_min_abscos > max_allowed_epip_abscos) {
				valid = false;
				break;
			}

			// add current image correspondences
			res.push_back(cur_epipolar_intersections);
		} else {
			// element corresponding to starting image in the defined order, just put a vector containing the starting intersection point coordinates
			std::vector<PolyLineGraph2D::plg_point> original_intersection;
			original_intersection.push_back(intersection_point_on_starting_image_plgp);
			res.push_back(original_intersection);
		}
	}
}

// Sets min_epip_angle to the minimum angle betweep epipolar and segment for all intersected segments
std::vector<PolyLineGraph2D::plg_point> PLGEdgeManagerClosestOnly::detect_epipolar_intersections_on_image_plgp_exclude_parallel_epipsegments(const int img_to_find_intersections_on, const glm::vec2 &starting_point_on_current_image, const float max_correspondence_detection_radius, const glm::vec3 &epipolar_line, float &max_epip_abscos) {
	std::vector<PolyLineGraph2D::plg_point> intersections;
	glm::vec2 cur_intersection;
	bool found;
	const vector<glm::vec4> &cur_segments = all_segments[img_to_find_intersections_on];
	float detection_dist_sq = pow(max_correspondence_detection_radius,2);

	max_epip_abscos = 0.0;
	float cur_epip_abscos;

	for(ulong i=0; i < cur_segments.size(); i++) {
		const glm::vec4 &segm = cur_segments[i];

#ifdef SELECT_CLOSE_POINTS_ENABLED
		single_intersect_segment_line_with_close_points_on_segm_in_detection_range(segm,epipolar_line, starting_point_on_current_image,detection_dist_sq, MAX_CLOSE_POINT_DISTANCESQ, cur_intersection, found);
#else
		single_intersect_segment_line_in_detection_range(segm,epipolar_line, starting_point_on_current_image,detection_dist_sq, cur_intersection,found);
#endif

		if (found) {
			// Compute angle
			cur_epip_abscos = abs(compute_anglecos(segm,epipolar_line));
			max_epip_abscos = cur_epip_abscos > max_epip_abscos ? cur_epip_abscos : max_epip_abscos;

			intersections.push_back(PolyLineGraph2D::plg_point(all_segments_polyline_ids[img_to_find_intersections_on][i],cur_intersection));
		}
	}

	return intersections;
}

std::vector<pair<PolyLineGraph2D::plg_point,std::vector<std::vector<PolyLineGraph2D::plg_point>>>> PLGEdgeManagerClosestOnly::detect_epipolar_intersections_for_all_starting_intersections_plgp_exclude_parallel_epipsegments(const int starting_image_id, const int starting_point_id, vector<PolyLineGraph2D::plg_point> &intersection_points_on_starting_image, const float correspondence_detection_range_multiplication_factor, const float min_correspondence_detection_radius, const float max_allowed_epip_abscos) {
	std::vector<pair<PolyLineGraph2D::plg_point,std::vector<std::vector<PolyLineGraph2D::plg_point>>>> res;
	bool valid;
	//std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>> result;

	bool tmp;
	glm::vec2 init_coords;
	get_2d_coordinates_of_point_on_image(sfmd,starting_image_id,starting_point_id,init_coords,tmp);

	for (const auto &intersection_point_on_starting_image : intersection_points_on_starting_image) {
		float max_correspondence_detection_radius = max(compute_2d_distance(init_coords,intersection_point_on_starting_image.plp.coords) * correspondence_detection_range_multiplication_factor, min_correspondence_detection_radius);
		std::vector<std::vector<PolyLineGraph2D::plg_point>> cur_res;
		detect_epipolar_intersections_plgp_exclude_parallel_epipsegments(starting_image_id, starting_point_id, intersection_point_on_starting_image, max_correspondence_detection_radius,max_allowed_epip_abscos,cur_res,valid);
		if(valid)
			res.push_back(make_pair(intersection_point_on_starting_image,cur_res));
	}

	return res;
}

std::vector<pair<PolyLineGraph2D::plg_point,std::vector<std::vector<PolyLineGraph2D::plg_point>>>> PLGEdgeManagerClosestOnly::detect_nearby_intersections_and_correspondences_plgp_exclude_parallel_epipsegments(const int starting_image_id, const int starting_point_id, const float starting_detection_dist, const float correspondence_detection_range_multiplication_factor, const float min_correspondence_detection_radius, const float max_allowed_epip_abscos) {
	std::vector<PolyLineGraph2D::plg_point> nearby_intersections = detect_nearby_edge_intersections_plgp(starting_image_id, starting_point_id, starting_detection_dist);
	std::vector<pair<PolyLineGraph2D::plg_point,std::vector<std::vector<PolyLineGraph2D::plg_point>>>> all_correspondences = detect_epipolar_intersections_for_all_starting_intersections_plgp_exclude_parallel_epipsegments(starting_image_id, starting_point_id, nearby_intersections, correspondence_detection_range_multiplication_factor,min_correspondence_detection_radius, max_allowed_epip_abscos);
	return all_correspondences;
}

std::vector<glm::vec2> PLGEdgeManagerClosestOnly::detect_nearby_edge_intersections(const int imgId, const int startingPoint_id, const float starting_detection_dist) {
	throw new NotImplementedException();
}
std::vector<std::vector<glm::vec2>> PLGEdgeManagerClosestOnly::detect_epipolar_intersections(const int starting_image_id, const int starting_point_id, const glm::vec2 &intersection_point_on_starting_image, const float max_correspondence_detection_radius) {
	throw new NotImplementedException();
}

std::vector<std::vector<std::vector<glm::vec2>>> PLGEdgeManagerClosestOnly::detect_epipolar_intersections_for_all_starting_intersections(const int starting_image_id, const int starting_point_id, const vector<glm::vec2> &intersection_points_on_starting_image, const float max_correspondence_detection_radius) {
	throw new NotImplementedException();
}

std::pair<std::vector<glm::vec2>,std::vector<std::vector<std::vector<glm::vec2>>>> PLGEdgeManagerClosestOnly::detect_nearby_intersections_and_correspondences(const int starting_image_id, const int starting_point_id, const float starting_detection_dist, const float max_correspondence_detection_radius) {
	throw new NotImplementedException();
}

std::vector<glm::vec2> PLGEdgeManagerClosestOnly::detect_epipolar_intersections_on_image(const int img_to_find_intersections_on, const glm::vec2 &starting_point_on_current_image, const float max_correspondence_detection_radius, const glm::vec3 &epipolar_line) {
	throw new NotImplementedException();
}
