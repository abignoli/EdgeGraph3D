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
#include <limits>
#include <set>

#include <glm.hpp>
#include "plg_edge_manager.hpp"

#include "SfMData.h"
#include "drawing_utilities.hpp"
#include "geometric_utilities.hpp"
#include "NotImplementedException.hpp"

PLGEdgeManager::PLGEdgeManager(const vector<Mat> &imgs, const SfMData &input_sfmd, const cv::Mat** input_all_fundamental_matrices, vector<PolyLineGraph2DHMapImpl> &input_plgs, const float starting_detection_dist, const float correspondence_detection_range_multiplication_factor) : EdgeManager(imgs,input_sfmd,input_all_fundamental_matrices), plgs(input_plgs), starting_detection_dist(starting_detection_dist), correspondence_detection_range_multiplication_factor(correspondence_detection_range_multiplication_factor), correspondence_detection_dist(starting_detection_dist * correspondence_detection_range_multiplication_factor), starting_detection_distsq(starting_detection_dist * starting_detection_dist), correspondence_detection_distsq(correspondence_detection_dist * correspondence_detection_dist) {
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

	for(const auto &plg: plgs)
		correspondence_plmaps.push_back(PolyLine2DMapSearch(plg,imgs[0].size(),correspondence_detection_dist));
}

void PLGEdgeManager::detect_edges(const vector<Mat> &imgs) {}

std::vector<std::vector<PolyLineGraph2D::plg_point>> PLGEdgeManager::detect_epipolar_intersections_plgp(const int starting_image_id, const int starting_point_id, const PolyLineGraph2D::plg_point &intersection_point_on_starting_image_plgp, const float max_correspondence_detection_radius) {
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

std::vector<PolyLineGraph2D::plg_point> PLGEdgeManager::detect_nearby_edge_intersections_plgp(const int imgId, const int startingPoint_id) {
	// get only closest polyline if close enough

	glm::vec2 startingPoint;
	bool tmp;

	// get 2D coordinates of starting point on starting image
	get_2d_coordinates_of_point_on_image(sfmd, imgId, startingPoint_id, startingPoint, tmp);

	const PolyLine2DMapSearch &plmap = correspondence_plmaps[imgId];

	return plmap.find_polylines_within_smaller_search_dist_with_plgps(startingPoint, starting_detection_dist);
}

std::vector<PolyLineGraph2D::plg_point> PLGEdgeManager::detect_epipolar_intersections_on_image_plgp(const int img_to_find_intersections_on, const glm::vec2 &starting_point_on_current_image, const float max_correspondence_detection_radius, const glm::vec3 &epipolar_line) {
	std::vector<PolyLineGraph2D::plg_point> intersections;

	float detection_dist_sq = max_correspondence_detection_radius * max_correspondence_detection_radius;

	const PolyLine2DMapSearch &plmap = correspondence_plmaps[img_to_find_intersections_on];
	set<ulong> potentially_valid_polylines = plmap.find_polylines_potentially_within_search_dist(starting_point_on_current_image);
	for(const auto pl_id : potentially_valid_polylines) {
		std::vector<PolyLineGraph2D::polyline::pl_point> pl_intersections = plgs[img_to_find_intersections_on].polylines[pl_id].intersect_line(epipolar_line);
		for(auto &plp: pl_intersections)
			if(squared_2d_distance(starting_point_on_current_image,plp.coords) <= detection_dist_sq)
				intersections.push_back(PolyLineGraph2D::plg_point(pl_id,plp));
	}

	return intersections;
}

vector<Mat> PLGEdgeManager::get_edge_images() {
	vector<Mat> edge_images;

	for(int i=0; i< all_segments.size(); i++) {
		Mat edge_image = draw_MultiColorComponents_PolyLineGraph_simplified(original_imgs[i], plgs[i]);
		edge_images.push_back(edge_image);
	}
	return edge_images;
}

vector<Mat> PLGEdgeManager::get_edge_images_original_background() {
	vector<Mat> edge_images;

	for(int i=0; i< all_segments.size(); i++) {
		Mat edge_image = Mat(original_imgs[0].rows,original_imgs[0].cols,CV_8UC3);
		draw_overlay_MultiColorComponents_PolyLineGraph_simplified(edge_image, plgs[i]);
		edge_images.push_back(edge_image);
	}
	return edge_images;
}

std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>> PLGEdgeManager::detect_epipolar_intersections_for_all_starting_intersections_plgp(const int starting_image_id, const int starting_point_id, const vector<PolyLineGraph2D::plg_point> &intersection_points_on_starting_image) {
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

std::pair<std::vector<PolyLineGraph2D::plg_point>,std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>>> PLGEdgeManager::detect_nearby_intersections_and_correspondences_plgp(const int starting_image_id, const int starting_point_id) {
	std::vector<PolyLineGraph2D::plg_point> nearby_intersections = detect_nearby_edge_intersections_plgp(starting_image_id, starting_point_id);
	std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>> all_correspondences = detect_epipolar_intersections_for_all_starting_intersections_plgp(starting_image_id, starting_point_id, nearby_intersections);
	return make_pair(nearby_intersections,all_correspondences);
}


std::vector<PolyLineGraph2D::plg_point> PLGEdgeManager::detect_epipolar_intersections_on_image_plgp(const int img_to_find_intersections_on, const glm::vec2 &starting_point_on_current_image, const float max_correspondence_detection_radius, const glm::vec3 &epipolar_line, const vector<ulong> &potentially_correspondent_polylines) {
	std::vector<PolyLineGraph2D::plg_point> intersections;

	float detection_dist_sq = max_correspondence_detection_radius*max_correspondence_detection_radius;

	const PolyLine2DMapSearch &plmap = correspondence_plmaps[img_to_find_intersections_on];
	for(const auto pl_id : potentially_correspondent_polylines) {
		std::vector<PolyLineGraph2D::polyline::pl_point> pl_intersections = plgs[img_to_find_intersections_on].polylines[pl_id].intersect_line(epipolar_line);
		for(auto &plp: pl_intersections)
			if(squared_2d_distance(starting_point_on_current_image,plp.coords) <= detection_dist_sq)
				intersections.push_back(PolyLineGraph2D::plg_point(pl_id,plp));
	}

	return intersections;
}


std::vector<std::vector<PolyLineGraph2D::plg_point>> PLGEdgeManager::detect_epipolar_intersections_plgp(const int starting_image_id, const int starting_point_id, const PolyLineGraph2D::plg_point &intersection_point_on_starting_image_plgp, const float max_correspondence_detection_radius, const vector<vector<ulong>> &potentially_correspondent_polylines) {
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
			cur_epipolar_intersections = detect_epipolar_intersections_on_image_plgp(cur_img, starting_point_on_cur_img, max_correspondence_detection_radius, epipolar_line, potentially_correspondent_polylines[i]);

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


std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>> PLGEdgeManager::detect_epipolar_intersections_for_all_starting_intersections_plgp(const int starting_image_id, const int starting_point_id, const vector<PolyLineGraph2D::plg_point> &intersection_points_on_starting_image, const vector<vector<ulong>> &potentially_correspondent_polylines) {
	std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>> result;

	bool tmp;
	glm::vec2 init_coords;
	get_2d_coordinates_of_point_on_image(sfmd,starting_image_id,starting_point_id,init_coords,tmp);

	for (const auto &intersection_point_on_starting_image : intersection_points_on_starting_image) {
		float max_correspondence_detection_radius = compute_2d_distance(init_coords,intersection_point_on_starting_image.plp.coords) * correspondence_detection_range_multiplication_factor;
		result.push_back(detect_epipolar_intersections_plgp(starting_image_id, starting_point_id, intersection_point_on_starting_image, max_correspondence_detection_radius, potentially_correspondent_polylines));
	}

	return result;
}

std::vector<std::pair<std::vector<PolyLineGraph2D::plg_point>,std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>>>> PLGEdgeManager::detect_nearby_intersections_and_correspondences_plgp(const int starting_point_id) {
	vector<vector<ulong>> potentially_correspondent_polylines;
	vector<vector<PolyLineGraph2D::plg_point>> starting_nearby_intersections;
	glm::vec2 startingPoint;
	bool tmp;

	ulong closest_segm;
	glm::vec2 projection;

	float cur_distsq;

	for(const auto starting_img_id : sfmd.camViewingPointN_[starting_point_id]) {
		// get 2D coordinates of starting point on starting image
		get_2d_coordinates_of_point_on_image(sfmd, starting_img_id, starting_point_id, startingPoint, tmp);
		set<ulong> cur_pcps_tmp = correspondence_plmaps[starting_img_id].find_polylines_potentially_within_search_dist(startingPoint);
		vector<ulong> cur_pcps;
		vector<PolyLineGraph2D::plg_point> cur_sni;
		for(const auto pl_id :cur_pcps_tmp) {
			cur_distsq = plgs[starting_img_id].polylines[pl_id].compute_distancesq(startingPoint, closest_segm, projection);
			if(cur_distsq <= starting_detection_distsq) {
				cur_pcps.push_back(pl_id);
				cur_sni.push_back(PolyLineGraph2D::plg_point(pl_id, closest_segm, projection));
			} else if(cur_distsq <= correspondence_detection_distsq)
				cur_pcps.push_back(pl_id);
		}
		potentially_correspondent_polylines.push_back(cur_pcps);
		starting_nearby_intersections.push_back(cur_sni);
	}

	std::vector<std::pair<std::vector<PolyLineGraph2D::plg_point>,std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>>>> res;
	for(int i=0; i < sfmd.camViewingPointN_[starting_point_id].size(); i++) {
		const auto starting_img_id = sfmd.camViewingPointN_[starting_point_id][i];
		std::vector<PolyLineGraph2D::plg_point> &nearby_intersections = starting_nearby_intersections[i];

		std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>> all_correspondences = detect_epipolar_intersections_for_all_starting_intersections_plgp(starting_img_id, starting_point_id, nearby_intersections, potentially_correspondent_polylines);
		res.push_back(make_pair(nearby_intersections,all_correspondences));
	}

	return res;
}

std::vector<std::pair<std::vector<PolyLineGraph2D::plg_point>,std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>>>> PLGEdgeManager::detect_nearby_intersections_and_correspondences_plgp_closest_starting_plgp(const int starting_point_id) {
	vector<vector<ulong>> potentially_correspondent_polylines;
	vector<vector<PolyLineGraph2D::plg_point>> starting_nearby_intersections;
	glm::vec2 startingPoint;
	bool tmp;

	ulong closest_segm;
	glm::vec2 projection;

	float cur_distsq;

	for(const auto starting_img_id : sfmd.camViewingPointN_[starting_point_id]) {
		bool found = false;
		float min_distsq = std::numeric_limits<float>::max();
		ulong min_pl_id;
		ulong min_closest_segm;
		glm::vec2 min_projection;

		// get 2D coordinates of starting point on starting image
		get_2d_coordinates_of_point_on_image(sfmd, starting_img_id, starting_point_id, startingPoint, tmp);
		set<ulong> cur_pcps_tmp = correspondence_plmaps[starting_img_id].find_polylines_potentially_within_search_dist(startingPoint);
		vector<ulong> cur_pcps;
		vector<PolyLineGraph2D::plg_point> cur_sni;
		for(const auto pl_id :cur_pcps_tmp) {
			cur_distsq = plgs[starting_img_id].polylines[pl_id].compute_distancesq(startingPoint, closest_segm, projection);
			if(cur_distsq <= starting_detection_distsq) {
				cur_pcps.push_back(pl_id);

				if(cur_distsq < min_distsq) {
					found = true;
					min_distsq = cur_distsq;
					min_pl_id = pl_id;
					min_closest_segm = closest_segm;
					min_projection=projection;
				}

			} else if(cur_distsq <= correspondence_detection_distsq)
				cur_pcps.push_back(pl_id);
		}

		if(found)
			cur_sni.push_back(PolyLineGraph2D::plg_point(min_pl_id, min_closest_segm, min_projection));

		potentially_correspondent_polylines.push_back(cur_pcps);
		starting_nearby_intersections.push_back(cur_sni);
	}

	std::vector<std::pair<std::vector<PolyLineGraph2D::plg_point>,std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>>>> res;
	for(int i=0; i < sfmd.camViewingPointN_[starting_point_id].size(); i++) {
		const auto starting_img_id = sfmd.camViewingPointN_[starting_point_id][i];
		std::vector<PolyLineGraph2D::plg_point> &nearby_intersections = starting_nearby_intersections[i];

		std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>> all_correspondences = detect_epipolar_intersections_for_all_starting_intersections_plgp(starting_img_id, starting_point_id, nearby_intersections, potentially_correspondent_polylines);
		res.push_back(make_pair(nearby_intersections,all_correspondences));
	}

	return res;
}

std::vector<glm::vec2> PLGEdgeManager::detect_nearby_edge_intersections(const int imgId, const int startingPoint_id, const float starting_detection_dist) {
	throw new NotImplementedException();
}
std::vector<std::vector<glm::vec2>> PLGEdgeManager::detect_epipolar_intersections(const int starting_image_id, const int starting_point_id, const glm::vec2 &intersection_point_on_starting_image, const float max_correspondence_detection_radius) {
	throw new NotImplementedException();
}

std::vector<std::vector<std::vector<glm::vec2>>> PLGEdgeManager::detect_epipolar_intersections_for_all_starting_intersections(const int starting_image_id, const int starting_point_id, const vector<glm::vec2> &intersection_points_on_starting_image, const float max_correspondence_detection_radius) {
	throw new NotImplementedException();
}

std::pair<std::vector<glm::vec2>,std::vector<std::vector<std::vector<glm::vec2>>>> PLGEdgeManager::detect_nearby_intersections_and_correspondences(const int starting_image_id, const int starting_point_id, const float starting_detection_dist, const float max_correspondence_detection_radius) {
	throw new NotImplementedException();
}

std::vector<glm::vec2> PLGEdgeManager::detect_epipolar_intersections_on_image(const int img_to_find_intersections_on, const glm::vec2 &starting_point_on_current_image, const float max_correspondence_detection_radius, const glm::vec3 &epipolar_line) {
	throw new NotImplementedException();
}
