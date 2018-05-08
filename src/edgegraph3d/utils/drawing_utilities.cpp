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

#include "drawing_utilities.hpp"

#include <opencv2/core/cvdef.h>
#include <opencv2/core/cvstd.hpp>
#include <opencv2/core/mat.inl.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc/types_c.h>
#include <opencv2/imgproc.hpp>
#include <stdlib.h>     /* srand, rand */
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>

#include "SfMData.h"
#include "geometric_utilities.hpp"
#include "glm.hpp"
#include "plg_analysis_utilities.hpp"
#include "convert_edge_images_pixel_to_segment.hpp"

using namespace std;
using namespace cv;

void draw_point(Mat &img, const Point2f &p, const Scalar &color, const int draw_radius) {
	if(draw_radius == 1)
		img.at<Vec3b>((int) p.y,(int) p.x) = Vec3b(color[0],color[1],color[2]);
	else
		cv::circle(img,p,draw_radius,color,DRAW_FILLED_CIRCLE);
}

void draw_point_glm(Mat &img, const glm::vec2 &p, const Scalar &color, const int draw_radius) {
	draw_point(img, Point2f(p[0],p[1]), color, draw_radius);
}

void draw_point(Mat &img, const Point2f &p, const int draw_radius) {
	Scalar color = generate_random_color();
	if(draw_radius == 1)
		img.at<Vec3b>((int) p.y,(int) p.x) = Vec3b(color[0],color[1],color[2]);
	else
		cv::circle(img,p,draw_radius,color,DRAW_FILLED_CIRCLE);
}

void draw_point_glm(Mat &img, const glm::vec2 &p, const int draw_radius) {
	draw_point(img, Point2f(p[0],p[1]), draw_radius);
}

void draw_points_glm(Mat &img, const vector<glm::vec2> &ps, const int draw_radius) {
	for(const auto &p: ps)
		draw_point_glm(img, p, draw_radius);
}

void draw_points_glm(Mat &img, const vector<glm::vec2> &ps, const Scalar &color, const int draw_radius) {
	for(const auto &p: ps)
		draw_point_glm(img, p, color, draw_radius);
}

void draw_points_glm(Mat &img, const vector<glm::vec2> &ps, const vector<Scalar> &colors, const int draw_radius) {
	for(int i=0; i < ps.size();i++) {
		const glm::vec2 p = ps[i];
		const Scalar &color = colors[i];
		draw_point_glm(img, p, color, draw_radius);
	}
}

void draw_reference_point_glm(Mat &img, const glm::vec2 &p, const Scalar &color) {
	draw_point(img, Point2f(p[0],p[1]), color, DRAW_REFERENCE_POINT_RADIUS);
}

void draw_reference_points_glm(Mat &img, const vector<glm::vec2> &ps, const vector<Scalar> &colors) {
	draw_points_glm(img, ps, colors, DRAW_REFERENCE_POINT_RADIUS);
}

void draw_intersection_point_glm(Mat &img, const glm::vec2 &p, const Scalar &color) {
	draw_point(img, Point2f(p[0],p[1]), color, DRAW_INTERSECTION_POINT_RADIUS);
}

void draw_segment_on_img(Mat &img, const glm::vec4 &segm, const Scalar &color) {
	cv::line(img, Point2f(segm[0],segm[1]), Point2f(segm[2],segm[3]), color);
}

void draw_segment_on_img(Mat &img, const glm::vec4 &segm, const Scalar &color, float thickness) {
	cv::line(img, Point2f(segm[0],segm[1]), Point2f(segm[2],segm[3]), color, thickness);
}

void draw_line_glm(Mat& img, const glm::vec3 &line, const Scalar& color) {
	int height = img.rows;
	int width = img.cols;

	Point2f pt1(0,- line[2] / line[1]);
	Point2f pt2(width,- (line[2] + width * line[0]) / line[1]);

	// cout << "Drawing epipolar from " << pt1 << " to " << pt2 << endl;

	cv::line(img, pt1, pt2, color);
}

void draw_lines_glm(Mat& img, const vector<glm::vec3> &lines, const vector<Scalar> &colors) {
	for(int i=0;i < lines.size();i++)
		draw_line_glm(img,lines[i],colors[i]);
}

void draw_circle_glm(Mat &img, const glm::vec2 &center,const int radius, const Scalar &color) {
	cv::circle(img,Point2f(center[0],center[1]),radius,color,DRAW_EMPTY_CIRCLE);
}

void draw_single_point_process(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<glm::vec2> &starting_intersections, const int correspondence_detection_radius, const std::vector< std::vector<std::vector<glm::vec2>>> &detected_correspondent_intersections) {
	Scalar starting_refpoint_color = generate_random_color();
	vector<Scalar> starting_intersection_colors = generate_random_colors(starting_intersections.size());
	draw_single_point_process(imgs, sfmd, all_fundamental_matrices, starting_img, reference_point_id, initial_detection_radius, starting_intersections, correspondence_detection_radius, detected_correspondent_intersections, starting_refpoint_color ,starting_intersection_colors);
}

vector<glm::vec2> convert_vec_plgpoint_to_vec2(const std::vector<PolyLineGraph2D::plg_point> &starting_intersections_plgp) {
	vector<glm::vec2> res;
	for(const auto &si : starting_intersections_plgp)
		res.push_back(si.plp.coords);
	return res;
}

vector<vector<vector<glm::vec2>>> convert_vecvecvec_plgpoint_to_vecvecvec2(const std::vector< std::vector<std::vector<PolyLineGraph2D::plg_point>>> &detected_correspondent_intersections_plgp) {
	vector<vector<vector<glm::vec2>>> res;
	for(const auto &cor : detected_correspondent_intersections_plgp) {
		vector<vector<glm::vec2>> cur_res;
		for(const auto &ccor : cor)
			cur_res.push_back(convert_vec_plgpoint_to_vec2(ccor));
		res.push_back(cur_res);
	}
	return res;
}


void draw_single_point_process(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<PolyLineGraph2D::plg_point> &starting_intersections_plgp, const std::vector< std::vector<std::vector<PolyLineGraph2D::plg_point>>> &detected_correspondent_intersections_plgp) {
	vector<glm::vec2> starting_intersections = convert_vec_plgpoint_to_vec2(starting_intersections_plgp);
	vector<vector<vector<glm::vec2>>> detected_correspondent_intersections = convert_vecvecvec_plgpoint_to_vecvecvec2(detected_correspondent_intersections_plgp);

	Scalar starting_refpoint_color = generate_random_color();
	vector<Scalar> starting_intersection_colors = generate_random_colors(starting_intersections.size());

	glm::vec2 p_cur2dcoords;
	const std::vector<int> &cams = sfmd.camViewingPointN_[reference_point_id];
	const std::vector<glm::vec2> &coords_on_cams = sfmd.point2DoncamViewingPoint_[reference_point_id];
	int cur_cam;
	int cam_index,starting_intersection_index;
	const int starting_intersections_num = starting_intersections.size();
	glm::vec3 epipolar_line;

	glm::vec2 refpoint_starting_image;
	bool found;
	get_2d_coordinates_of_point_on_image(sfmd,starting_img,reference_point_id,refpoint_starting_image,found);


	// Draw reference point and, only on starting image , corresponding detection radius and starting intersections
	for(cam_index=0; cam_index < cams.size(); cam_index++) {
		cur_cam = cams[cam_index];
		p_cur2dcoords = coords_on_cams[cam_index];

		// Draw reference starting point on each image where it's visible
		draw_reference_point_glm(imgs[cur_cam], p_cur2dcoords, starting_refpoint_color);

		if(cur_cam == starting_img) {
			// Draw starting detection radius
			draw_circle_glm(imgs[starting_img], p_cur2dcoords,initial_detection_radius, starting_refpoint_color);

			for(starting_intersection_index=0; starting_intersection_index < starting_intersections_num; starting_intersection_index++) {
				// Draw starting intersection on first image
				draw_intersection_point_glm(imgs[starting_img], starting_intersections[starting_intersection_index], starting_intersection_colors[starting_intersection_index]);
			}
		} else {
			// Get Fundamental Matrix
			cv::Mat &F = all_fundamental_matrices[starting_img][cur_cam];

			// Draw epipolar line on refpoint
			if(computeCorrespondEpilineSinglePoint(refpoint_starting_image, F, epipolar_line,1))
				draw_line_glm(imgs[cur_cam], epipolar_line, starting_refpoint_color);

			for(starting_intersection_index=0; starting_intersection_index < starting_intersections_num; starting_intersection_index++) {
				//Draw epipolar line
				if(computeCorrespondEpilineSinglePoint(starting_intersections[starting_intersection_index], F, epipolar_line,1))
					draw_line_glm(imgs[cur_cam], epipolar_line, starting_intersection_colors[starting_intersection_index]);

				for(const auto &correspondence_point : detected_correspondent_intersections[starting_intersection_index][cam_index]) {
					// Draw starting intersection on first image
					draw_intersection_point_glm(imgs[cur_cam], correspondence_point, starting_intersection_colors[starting_intersection_index]);
				}
			}
		}
	}}

void draw_single_point_process(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<pair<PolyLineGraph2D::plg_point,std::vector<std::vector<PolyLineGraph2D::plg_point>>>> &nearby_matches) {

	Scalar starting_refpoint_color = generate_random_color();
	vector<Scalar> starting_intersection_colors = generate_random_colors(nearby_matches.size());

	glm::vec2 p_cur2dcoords;
	const std::vector<int> &cams = sfmd.camViewingPointN_[reference_point_id];
	const std::vector<glm::vec2> &coords_on_cams = sfmd.point2DoncamViewingPoint_[reference_point_id];
	int cur_cam;
	int cam_index,starting_intersection_index;
	const int starting_intersections_num = nearby_matches.size();
	glm::vec3 epipolar_line;

	glm::vec2 refpoint_starting_image;
	bool found;
	get_2d_coordinates_of_point_on_image(sfmd,starting_img,reference_point_id,refpoint_starting_image,found);


	// Draw reference point and, only on starting image , corresponding detection radius and starting intersections
	for(cam_index=0; cam_index < cams.size(); cam_index++) {
		cur_cam = cams[cam_index];
		p_cur2dcoords = coords_on_cams[cam_index];

		// Draw reference starting point on each image where it's visible
		draw_reference_point_glm(imgs[cur_cam], p_cur2dcoords, starting_refpoint_color);

		if(cur_cam == starting_img) {
			// Draw starting detection radius
			draw_circle_glm(imgs[starting_img], p_cur2dcoords,initial_detection_radius, starting_refpoint_color);

			for(starting_intersection_index=0; starting_intersection_index < starting_intersections_num; starting_intersection_index++) {
				// Draw starting intersection on first image
				draw_intersection_point_glm(imgs[starting_img], nearby_matches[starting_intersection_index].first.plp.coords, starting_intersection_colors[starting_intersection_index]);
			}
		} else {
			// Get Fundamental Matrix
			cv::Mat &F = all_fundamental_matrices[starting_img][cur_cam];

			for(starting_intersection_index=0; starting_intersection_index < starting_intersections_num; starting_intersection_index++) {
				//Draw epipolar line
				if(computeCorrespondEpilineSinglePoint(nearby_matches[starting_intersection_index].first.plp.coords, F, epipolar_line,1))
					draw_line_glm(imgs[cur_cam], epipolar_line, starting_intersection_colors[starting_intersection_index]);

				for(const auto &correspondence_point : nearby_matches[starting_intersection_index].second[cam_index]) {
					// Draw starting intersection on first image
					draw_intersection_point_glm(imgs[cur_cam], correspondence_point.plp.coords, starting_intersection_colors[starting_intersection_index]);
				}
			}
		}
	}
}

void draw_single_point_process_no_epilines(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<pair<PolyLineGraph2D::plg_point,std::vector<std::vector<PolyLineGraph2D::plg_point>>>> &nearby_matches) {

	Scalar starting_refpoint_color = generate_random_color();
	vector<Scalar> starting_intersection_colors = generate_random_colors(nearby_matches.size());

	glm::vec2 p_cur2dcoords;
	const std::vector<int> &cams = sfmd.camViewingPointN_[reference_point_id];
	const std::vector<glm::vec2> &coords_on_cams = sfmd.point2DoncamViewingPoint_[reference_point_id];
	int cur_cam;
	int cam_index,starting_intersection_index;
	const int starting_intersections_num = nearby_matches.size();
	glm::vec3 epipolar_line;

	glm::vec2 refpoint_starting_image;
	bool found;
	get_2d_coordinates_of_point_on_image(sfmd,starting_img,reference_point_id,refpoint_starting_image,found);


	// Draw reference point and, only on starting image , corresponding detection radius and starting intersections
	for(cam_index=0; cam_index < cams.size(); cam_index++) {
		cur_cam = cams[cam_index];
		p_cur2dcoords = coords_on_cams[cam_index];

		if(cur_cam == starting_img) {

			for(starting_intersection_index=0; starting_intersection_index < starting_intersections_num; starting_intersection_index++) {
				// Draw starting intersection on first image
				draw_intersection_point_glm(imgs[starting_img], nearby_matches[starting_intersection_index].first.plp.coords, starting_intersection_colors[starting_intersection_index]);
			}
		} else {
			// Get Fundamental Matrix
			cv::Mat &F = all_fundamental_matrices[starting_img][cur_cam];

			// Draw epipolar line on refpoint
			computeCorrespondEpilineSinglePoint(refpoint_starting_image, F, epipolar_line,1);

			for(starting_intersection_index=0; starting_intersection_index < starting_intersections_num; starting_intersection_index++) {
				//Draw epipolar line
				computeCorrespondEpilineSinglePoint(nearby_matches[starting_intersection_index].first.plp.coords, F, epipolar_line,1);

				for(const auto &correspondence_point : nearby_matches[starting_intersection_index].second[cam_index]) {
					// Draw starting intersection on first image
					draw_intersection_point_glm(imgs[cur_cam], correspondence_point.plp.coords, starting_intersection_colors[starting_intersection_index]);
				}
			}
		}
	}
}

void draw_single_point_process_with_segment_output(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<SegmentEdgeManager::segment_point> &starting_intersections, const int correspondence_detection_radius, const std::vector< std::vector<std::vector<SegmentEdgeManager::segment_point>>> &detected_correspondent_intersections) {
	Scalar starting_refpoint_color = generate_random_color();
	vector<Scalar> starting_intersection_colors = generate_random_colors(starting_intersections.size());
	draw_single_point_process_with_segment_output(imgs, sfmd, all_fundamental_matrices, starting_img, reference_point_id, initial_detection_radius, starting_intersections, correspondence_detection_radius, detected_correspondent_intersections, starting_refpoint_color ,starting_intersection_colors);
}

/**
 * detected_correspondent_intersections : for each starting intersection point, for each image where reference point is visible except starting image, each correspondent intersection found
 */
void draw_single_point_process(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<glm::vec2> &starting_intersections, const int correspondence_detection_radius, const std::vector< std::vector<std::vector<glm::vec2>>> &detected_correspondent_intersections, const Scalar &starting_refpoint_color ,const vector<Scalar> &starting_intersection_colors) {
	glm::vec2 p_cur2dcoords;
	const std::vector<int> &cams = sfmd.camViewingPointN_[reference_point_id];
	const std::vector<glm::vec2> &coords_on_cams = sfmd.point2DoncamViewingPoint_[reference_point_id];
	int cur_cam;
	int cam_index,starting_intersection_index;
	const int starting_intersections_num = starting_intersections.size();
	glm::vec3 epipolar_line;

	glm::vec2 refpoint_starting_image;
	bool found;
	get_2d_coordinates_of_point_on_image(sfmd,starting_img,reference_point_id,refpoint_starting_image,found);


	// Draw reference point and, only on starting image , corresponding detection radius and starting intersections
	for(cam_index=0; cam_index < cams.size(); cam_index++) {
		cur_cam = cams[cam_index];
		p_cur2dcoords = coords_on_cams[cam_index];

		// Draw reference starting point on each image where it's visible
		draw_reference_point_glm(imgs[cur_cam], p_cur2dcoords, starting_refpoint_color);

		if(cur_cam == starting_img) {
			// Draw starting detection radius
			draw_circle_glm(imgs[starting_img], p_cur2dcoords,initial_detection_radius, starting_refpoint_color);

			for(starting_intersection_index=0; starting_intersection_index < starting_intersections_num; starting_intersection_index++) {
				// Draw starting intersection on first image
				draw_intersection_point_glm(imgs[starting_img], starting_intersections[starting_intersection_index], starting_intersection_colors[starting_intersection_index]);
			}
		} else {
			// Draw max range detection radius on other images
			draw_circle_glm(imgs[cur_cam], p_cur2dcoords,correspondence_detection_radius, starting_refpoint_color);

			// Get Fundamental Matrix
			cv::Mat &F = all_fundamental_matrices[starting_img][cur_cam];

			// Draw epipolar line on refpoint
			if(computeCorrespondEpilineSinglePoint(refpoint_starting_image, F, epipolar_line,1))
				draw_line_glm(imgs[cur_cam], epipolar_line, starting_refpoint_color);

			for(starting_intersection_index=0; starting_intersection_index < starting_intersections_num; starting_intersection_index++) {
				//Draw epipolar line
				if(computeCorrespondEpilineSinglePoint(starting_intersections[starting_intersection_index], F, epipolar_line,1))
					draw_line_glm(imgs[cur_cam], epipolar_line, starting_intersection_colors[starting_intersection_index]);

				for(const auto &correspondence_point : detected_correspondent_intersections[starting_intersection_index][cam_index]) {
					// Draw starting intersection on first image
					draw_intersection_point_glm(imgs[cur_cam], correspondence_point, starting_intersection_colors[starting_intersection_index]);
				}
			}
		}
	}
}

/**
 * detected_correspondent_intersections : for each starting intersection point, for each image where reference point is visible except starting image, each correspondent intersection found
 */
void draw_single_point_process_with_segment_output(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<SegmentEdgeManager::segment_point> &starting_intersections, const int correspondence_detection_radius, const std::vector< std::vector<std::vector<SegmentEdgeManager::segment_point>>> &detected_correspondent_intersections, const Scalar &starting_refpoint_color ,const vector<Scalar> &starting_intersection_colors) {
	glm::vec2 p_cur2dcoords;
	const std::vector<int> &cams = sfmd.camViewingPointN_[reference_point_id];
	const std::vector<glm::vec2> &coords_on_cams = sfmd.point2DoncamViewingPoint_[reference_point_id];
	int cur_cam;
	int cam_index,starting_intersection_index;
	const int starting_intersections_num = starting_intersections.size();
	glm::vec3 epipolar_line;

	glm::vec2 refpoint_starting_image;
	bool found;
	get_2d_coordinates_of_point_on_image(sfmd,starting_img,reference_point_id,refpoint_starting_image,found);


	// Draw reference point and, only on starting image , corresponding detection radius and starting intersections
	for(cam_index=0; cam_index < cams.size(); cam_index++) {
		cur_cam = cams[cam_index];
		p_cur2dcoords = coords_on_cams[cam_index];

		// Draw reference starting point on each image where it's visible
		draw_reference_point_glm(imgs[cur_cam], p_cur2dcoords, starting_refpoint_color);

		if(cur_cam == starting_img) {
			// Draw starting detection radius
			draw_circle_glm(imgs[starting_img], p_cur2dcoords,initial_detection_radius, starting_refpoint_color);

			for(starting_intersection_index=0; starting_intersection_index < starting_intersections_num; starting_intersection_index++) {
				// Draw starting intersection on first image
				draw_intersection_point_glm(imgs[starting_img], starting_intersections[starting_intersection_index].coords, starting_intersection_colors[starting_intersection_index]);
				draw_segment_on_img(imgs[starting_img], starting_intersections[starting_intersection_index].segment_extremes, starting_intersection_colors[starting_intersection_index]);
			}
		} else {
			// Draw max range detection radius on other images
			draw_circle_glm(imgs[cur_cam], p_cur2dcoords,correspondence_detection_radius, starting_refpoint_color);

			// Get Fundamental Matrix
			cv::Mat &F = all_fundamental_matrices[starting_img][cur_cam];

			// Draw epipolar line on refpoint
			if(computeCorrespondEpilineSinglePoint(refpoint_starting_image, F, epipolar_line,1))
				draw_line_glm(imgs[cur_cam], epipolar_line, starting_refpoint_color);

			for(starting_intersection_index=0; starting_intersection_index < starting_intersections_num; starting_intersection_index++) {
				//Draw epipolar line
				if(computeCorrespondEpilineSinglePoint(starting_intersections[starting_intersection_index].coords, F, epipolar_line,1))
					draw_line_glm(imgs[cur_cam], epipolar_line, starting_intersection_colors[starting_intersection_index]);

				for(const auto &correspondence_point : detected_correspondent_intersections[starting_intersection_index][cam_index]) {
					// Draw starting intersection on first image
					draw_intersection_point_glm(imgs[cur_cam], correspondence_point.coords, starting_intersection_colors[starting_intersection_index]);
					draw_segment_on_img(imgs[cur_cam], correspondence_point.segment_extremes, starting_intersection_colors[starting_intersection_index]);
				}
			}
		}
	}
}

void draw_single_point_process_and_consensus(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<SegmentEdgeManager::segment_point> &starting_intersections_w_segments, const int correspondence_detection_radius, const std::vector< std::vector<std::vector<SegmentEdgeManager::segment_point>>> &detected_correspondent_intersections_w_segments, const vector<std::tuple<glm::vec3,vector<glm::vec2>,vector<int>>> &found_3d_points, const Scalar &consensus_match_color) {
	std::vector<glm::vec2> starting_intersections = intersections_remove_segments(starting_intersections_w_segments);
	std::vector< std::vector<std::vector<glm::vec2>>> detected_correspondent_intersections = correspondences_remove_segments(detected_correspondent_intersections_w_segments);
	draw_single_point_process(imgs, sfmd, all_fundamental_matrices, starting_img, reference_point_id, initial_detection_radius, starting_intersections, correspondence_detection_radius, detected_correspondent_intersections);
	draw_consensus_matched_points(imgs,found_3d_points,consensus_match_color);
}

void draw_single_point_process_and_consensus(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<glm::vec2> &starting_intersections, const int correspondence_detection_radius, const std::vector< std::vector<std::vector<glm::vec2>>> &detected_correspondent_intersections, const vector<std::tuple<glm::vec3,vector<glm::vec2>,vector<int>>> &found_3d_points, const Scalar &consensus_match_color) {
	draw_single_point_process(imgs, sfmd, all_fundamental_matrices, starting_img, reference_point_id, initial_detection_radius, starting_intersections, correspondence_detection_radius, detected_correspondent_intersections);
	draw_consensus_matched_points(imgs,found_3d_points,consensus_match_color);
}

void draw_point_on_imgs(vector<Mat> &imgs,const vector<int> &cam_ids, const vector<glm::vec2> &coords) {
	Scalar col = generate_random_color();
	for(int i=0; i < cam_ids.size();i++)
		draw_intersection_point_glm(imgs[cam_ids[i]],coords[i],col);
}

void draw_new_consensus_points(vector<Mat> &imgs, const vector<std::tuple<glm::vec3,vector<glm::vec2>,vector<int>>> &found_3d_points) {
	for(const auto &p3d: found_3d_points)
		draw_point_on_imgs(imgs,std::get<2>(p3d),std::get<1>(p3d));
}

void draw_refpoint_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const int point_id, const Scalar &point_color) {
	const std::vector<int> &cams = sfmd.camViewingPointN_[point_id];
	const std::vector<glm::vec2> &coords_on_cams = sfmd.point2DoncamViewingPoint_[point_id];
	int cam_id;
	glm::vec2 p_cur2dcoords;

	for(int cam_index=0; cam_index<cams.size();cam_index++) {
		cam_id = cams[cam_index];
		p_cur2dcoords = coords_on_cams[cam_index];
		draw_reference_point_glm(imgs[cam_id], p_cur2dcoords, point_color);
	}
}

void draw_refpoints_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const vector<Scalar> &point_colors) {
	for(int point_id=0; point_id < sfmd.numPoints_; point_id++) {
		draw_refpoint_on_imgs(imgs,sfmd, point_id, point_colors[point_id]);
	}
}

void draw_refpoint_with_circle_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const int point_id, const Scalar &point_color, const float radius) {
	const std::vector<int> &cams = sfmd.camViewingPointN_[point_id];
	const std::vector<glm::vec2> &coords_on_cams = sfmd.point2DoncamViewingPoint_[point_id];
	int cam_id;
	glm::vec2 p_cur2dcoords;

	for(int cam_index=0; cam_index<cams.size();cam_index++) {
		cam_id = cams[cam_index];
		p_cur2dcoords = coords_on_cams[cam_index];
		draw_reference_point_glm(imgs[cam_id], p_cur2dcoords, point_color);
		draw_circle_glm(imgs[cam_id],p_cur2dcoords,radius,point_color);
	}
}

void draw_refpoints_with_circle_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const vector<Scalar> &point_colors, const float radius){
	for(int point_id=0; point_id < sfmd.numPoints_; point_id++) {
		draw_refpoint_with_circle_on_imgs(imgs,sfmd, point_id, point_colors[point_id],radius);
	}
}

void draw_refpoint_with_circles_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const int point_id, const Scalar &point_color, const float radius1, const float radius2) {
	const std::vector<int> &cams = sfmd.camViewingPointN_[point_id];
	const std::vector<glm::vec2> &coords_on_cams = sfmd.point2DoncamViewingPoint_[point_id];
	int cam_id;
	glm::vec2 p_cur2dcoords;

	for(int cam_index=0; cam_index<cams.size();cam_index++) {
		cam_id = cams[cam_index];
		p_cur2dcoords = coords_on_cams[cam_index];
		draw_reference_point_glm(imgs[cam_id], p_cur2dcoords, point_color);
		draw_circle_glm(imgs[cam_id],p_cur2dcoords,radius1,point_color);
		draw_circle_glm(imgs[cam_id],p_cur2dcoords,radius2,point_color);
	}
}

void draw_refpoints_with_circles_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const vector<Scalar> &point_colors, const float radius1, const float radius2) {
	for(int point_id=0; point_id < sfmd.numPoints_; point_id++) {
		draw_refpoint_with_circles_on_imgs(imgs,sfmd, point_id, point_colors[point_id],radius1,radius2);
}
}

void draw_setofrefpoints_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const vector<int> refpoints, const vector<Scalar> &point_colors) {
	for(const auto point_id : refpoints) {
		draw_refpoint_on_imgs(imgs,sfmd, point_id, point_colors[point_id]);
	}
}

void draw_refpoint_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const int point_id) {
	draw_refpoint_on_imgs(imgs,sfmd, point_id, generate_random_color());
}

void draw_refpoints_on_imgs(vector<Mat> &imgs,const SfMData &sfmd) {
	for(int point_id=0; point_id < sfmd.numPoints_; point_id++)
		draw_refpoint_on_imgs(imgs,sfmd, point_id);
}

void draw_point_projections(vector<Mat> &imgs, const vector<glm::vec2> &coords, const vector<int> &cameras, const Scalar &point_color)
{
	for(int i=0; i < coords.size(); i++)
		draw_reference_point_glm(imgs[cameras[i]], coords[i], point_color);
}

void draw_point_projections(vector<Mat> &imgs, const vector<glm::vec2> &coords, const vector<int> &cameras)
{
	const Scalar point_color = generate_random_color();
	for(int i=0; i < coords.size(); i++)
		draw_reference_point_glm(imgs[cameras[i]], coords[i], point_color);
}

void draw_3dpoint_on_imgs(vector<Mat> &imgs, const std::tuple<glm::vec3, vector<glm::vec2>, vector<int>> &p3d, const Scalar &point_color) {
	draw_point_projections(imgs,get<1>(p3d),get<2>(p3d),point_color);
}

void draw_3dpoint_on_imgs(vector<Mat> &imgs, const std::tuple<glm::vec3, vector<glm::vec2>, vector<int>> &p3d) {
	draw_3dpoint_on_imgs(imgs,p3d,generate_random_color());
}

void draw_3dpoints_on_imgs(vector<Mat> &imgs, const vector<std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>> &p3ds) {
	for(const auto &p3d: p3ds)
		draw_3dpoint_on_imgs(imgs,p3d,generate_random_color());
}

void draw_point_projections(vector<Mat> &imgs, const vector<PolyLineGraph2D::plg_point> &coords, const vector<int> &cameras, const Scalar &point_color)
{
	for(int i=0; i < coords.size(); i++)
		draw_reference_point_glm(imgs[cameras[i]], coords[i].plp.coords, point_color);
}

void draw_point_projections(vector<Mat> &imgs, const vector<PolyLineGraph2D::plg_point> &coords, const vector<int> &cameras)
{
	const Scalar point_color = generate_random_color();
	for(int i=0; i < coords.size(); i++)
		draw_reference_point_glm(imgs[cameras[i]], coords[i].plp.coords, point_color);
}

void draw_3dpoint_on_imgs(vector<Mat> &imgs, const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &p3d, const Scalar &point_color) {
	draw_point_projections(imgs,get<1>(p3d),get<2>(p3d),point_color);
}

void draw_3dpoint_on_imgs(vector<Mat> &imgs, const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &p3d) {
	draw_3dpoint_on_imgs(imgs,p3d,generate_random_color());
}

void draw_3dpoints_on_imgs(vector<Mat> &imgs, const vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &p3ds) {
	for(const auto &p3d: p3ds)
		draw_3dpoint_on_imgs(imgs,p3d,generate_random_color());
}

void draw_refpoint_with_circle_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const int point_id, const float radius) {
	draw_refpoint_with_circle_on_imgs(imgs,sfmd, point_id, generate_random_color(), radius);
}

void draw_refpoints_with_circle_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const float radius) {
	for(int point_id=0; point_id < sfmd.numPoints_; point_id++)
		draw_refpoint_with_circle_on_imgs(imgs,sfmd, point_id, radius);
}

void draw_refpoint_with_circles_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const int point_id, const float radius1, const float radius2) {
	draw_refpoint_with_circles_on_imgs(imgs, sfmd, point_id, generate_random_color(), radius1, radius2);
}

void draw_refpoints_with_circles_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const float radius1, const float radius2) {
	for(int point_id=0; point_id < sfmd.numPoints_; point_id++)
		draw_refpoint_with_circles_on_imgs(imgs,sfmd, point_id, radius1,radius2);
}

void draw_point_epipolars_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, glm::vec2 p, const int starting_img, const Scalar &point_color) {
	draw_reference_point_glm(imgs[starting_img], p, point_color);
	glm::vec3 epipolar_line;

	for(int cam_index=0; cam_index<sfmd.numCameras_;cam_index++)
		if(cam_index!=starting_img) {
			cv::Mat &F = all_fundamental_matrices[starting_img][cam_index];
			if(computeCorrespondEpilineSinglePoint(p, F, epipolar_line, 1))
				draw_line_glm(imgs[cam_index], epipolar_line, point_color);
		}
}

void draw_refpoint_epipolars_on_imgs(vector<Mat> &imgs,const SfMData &sfmd,const cv::Mat** all_fundamental_matrices, const int point_id, const int starting_img, const Scalar &point_color) {
	draw_refpoint_on_imgs(imgs,sfmd, point_id, point_color);
	glm::vec2 point_coords;
	bool found;
	get_2d_coordinates_of_point_on_image(sfmd, starting_img, point_id, point_coords, found);
	draw_point_epipolars_on_imgs(imgs,sfmd, all_fundamental_matrices, point_coords, starting_img, point_color);
}

void draw_setofrefpoints_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const vector<int> refpoints) {
	vector<Scalar> colors = generate_random_colors(sfmd.numPoints_);
	draw_setofrefpoints_on_imgs(imgs,sfmd,refpoints,colors);
}

Scalar generate_random_color() {
	return Scalar(rand() % 256,rand() % 256,rand() % 256);
}

vector<Scalar> generate_random_colors(int size) {
	vector<Scalar> result;
	for(int i=0;i<size;i++)
		result.push_back(generate_random_color());
	return result;
}

void draw_img_pair_refpoints(vector<Mat> &imgs,const SfMData &sfmd,const int i,const int j) {
	vector<set<int>> points_on_images = get_point_sets_on_images(sfmd);
	set<int> points_on_both = find_points_on_both_images(points_on_images, i, j);
	vector<int> points_on_both_vec;
	std::copy(points_on_both.begin(), points_on_both.end(), std::back_inserter(points_on_both_vec));
	vector<glm::vec2> points_i = get_vector_of_2d_coordinates_of_point_on_image(sfmd, i, points_on_both_vec);
	vector<glm::vec2> points_j = get_vector_of_2d_coordinates_of_point_on_image(sfmd, j, points_on_both_vec);

	const vector<Scalar> colors = generate_random_colors(points_on_both.size());

	draw_reference_points_glm(imgs[i], points_i, colors);
	draw_reference_points_glm(imgs[j], points_j, colors);
}

void draw_img_pair_epipolars_refpoints(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices,const int i,const int j) {
	vector<set<int>> points_on_images = get_point_sets_on_images(sfmd);
	set<int> points_on_both = find_points_on_both_images(points_on_images, i, j);
	vector<int> points_on_both_vec;
	std::copy(points_on_both.begin(), points_on_both.end(), std::back_inserter(points_on_both_vec));
	vector<glm::vec2> points_i = get_vector_of_2d_coordinates_of_point_on_image(sfmd, i, points_on_both_vec);
	vector<glm::vec2> points_j = get_vector_of_2d_coordinates_of_point_on_image(sfmd, j, points_on_both_vec);


	const Mat &F1 = all_fundamental_matrices[j][i];
	vector<glm::vec3> lines_i = computeCorrespondEpilineMultiplePoints(points_j, F1, 1);

	const Mat &F2 = all_fundamental_matrices[i][j];
	vector<glm::vec3> lines_j = computeCorrespondEpilineMultiplePoints(points_i, F2, 1);

	const vector<Scalar> colors = generate_random_colors(points_on_both.size());

	draw_reference_points_glm(imgs[i], points_i, colors);
	draw_reference_points_glm(imgs[j], points_j, colors);

	draw_lines_glm(imgs[i], lines_i, colors);
	draw_lines_glm(imgs[j], lines_j, colors);
}

/**
 * Images are arranged in a squared display, with no scaling.
 */
Mat getUnifiedSquareImage(const vector<Mat>& imgs, const int image_type) {
    float nImgs=imgs.size();
    int   rows=ceil(sqrt(nImgs));     // You can set this explicitly
    int   cols=ceil(nImgs/rows); // You can set this explicitly
    Size totalSize(cols*imgs[0].cols,rows*imgs[0].rows);
    return getUnifiedSquareImage(imgs, totalSize, image_type);
}

/**
 * Images are arranged in a squared display, of the specified size.
 */
Mat getUnifiedSquareImage(const vector<Mat>& imgs,const Size &totalSize, const int image_type) {
	Mat scaledImg;
    float nImgs=imgs.size();
    int   rows=ceil(sqrt(nImgs));     // You can set this explicitly
    int   cols=ceil(nImgs/rows); // You can set this explicitly
    Size cellSize(totalSize.width / cols, totalSize.height / rows);

    Size scaledSize;
    double fx,fy,scale;

    int resultImgW=cellSize.width*cols;
    int resultImgH=cellSize.height*rows;

    Mat resultImg=Mat::zeros(resultImgH,resultImgW,image_type);
    int ind=0;

    for(int i=0;i<rows;i++)
    {
        for(int j=0;j<cols;j++)
        {
            if(ind<imgs.size())
            {
				fx = float(cellSize.width) / imgs[ind].cols;
				fy = float(cellSize.height) / imgs[ind].rows;
				scale = fx < fy ? fx : fy;
				scaledSize.width = imgs[ind].cols * scale;
				scaledSize.height = imgs[ind].rows * scale;

				resize(imgs[ind], scaledImg, scaledSize, 0, 0, INTER_LINEAR);


				int cell_row=i*cellSize.height;
				int cell_col=j*cellSize.width;

				scaledImg.copyTo(resultImg(Range(cell_row,cell_row+scaledImg.rows),Range(cell_col,cell_col+scaledImg.cols)));
            }
            ind++;
        }
    }

    return resultImg;
}

// --------------------------------------------------------------
// Function to draw several images to one image.
// Small images draws into cells of size cellSize.
// If image larger than size of cell ot will be trimmed.
// If image smaller than cellSize there will be gap between cells.
// --------------------------------------------------------------
char showImages(string title,const vector<Mat>& imgs,const Size &totalSize, const int image_type)
{
	char k=0;

	Mat resultImg = getUnifiedSquareImage(imgs,totalSize,image_type);

    imshow(title,resultImg);
    //while(true)
    	k=waitKey(0);
    return k;
}

void draw_consensus_matched_point_on_img(const Mat &img,const glm::vec2 &coords, const Scalar &color) {
	cv::circle(img,Point2f(coords[0],coords[1]),DRAW_NEW_MATCHED_POINT_RADIUS,color,DRAW_EMPTY_CIRCLE);
}

void draw_consensus_matched_point(const vector<Mat>& imgs,const vector<int> &cam_ids, const vector<glm::vec2> &cam_coords, const Scalar &color) {
	for(int i=0; i< cam_ids.size();i++) {
		draw_consensus_matched_point_on_img(imgs[cam_ids[i]], cam_coords[i], color);
	}
}

void draw_consensus_matched_point(const vector<Mat>& imgs,const vector<int> &cam_ids, const vector<glm::vec2> &cam_coords) {
	draw_consensus_matched_point(imgs,cam_ids, cam_coords, generate_random_color());
}

void draw_consensus_matched_points(const vector<Mat>& imgs,const vector<std::tuple<glm::vec3,vector<glm::vec2>,vector<int>>> &new_points) {
	for(const auto &new_point : new_points) {
		draw_consensus_matched_point(imgs,std::get<2>(new_point),std::get<1>(new_point));
	}
}

void draw_consensus_matched_points(const vector<Mat>& imgs,const vector<std::tuple<glm::vec3,vector<glm::vec2>,vector<int>>> &new_points, const Scalar &color) {
	for(const auto &new_point : new_points) {
		draw_consensus_matched_point(imgs,std::get<2>(new_point),std::get<1>(new_point),color);
	}
}

void draw_segments_on_image(Mat &img, const vector<glm::vec4> &segments) {
	draw_segments_on_image(img, convert_vec4_Vec4f(segments));
}

Mat draw_segments_on_newimage(const Size &sz, const vector<glm::vec4> &segments, const Scalar &colorbg, const Scalar &colorlines) {
	int rows=sz.height;
	int cols=sz.width;
	Mat out = Mat(rows,cols,CV_8UC3);
	out.setTo(colorbg);
	for(const auto &s:segments)
		cv::line(out, cv::Point2f(s[0],s[1]), cv::Point2f(s[2],s[3]), colorlines);

	return out;
}

Mat draw_segments_on_newimage_with_extremes(const Size &sz, const vector<glm::vec4> &segments, const Scalar &colorbg, const Scalar &colorlines, const Scalar &colorstart, const Scalar &colorend) {
	int rows=sz.height;
	int cols=sz.width;
	Mat out = Mat(rows,cols,CV_8UC3);
	out.setTo(colorbg);
	for(const auto &s:segments) {
		if(s[0] >= 0 && s[0] < cols && s[1] >= 0 &&  s[1] < rows)
			out.at<Vec3b>(Point(int(s[0]),int(s[1])))=Vec3b(colorstart[0],colorstart[1],colorstart[2]);
		if(s[2] >= 0 && s[2] < cols && s[3] >= 0 && s[3] < rows)
			out.at<Vec3b>(Point(int(s[2]),int(s[3])))=Vec3b(colorend[0],colorend[1],colorend[2]);
		cv::line(out, cv::Point2f(s[0],s[1]), cv::Point2f(s[2],s[3]), colorlines);
	}
	return out;
}

void draw_segments_on_image_rnd_color(Mat &img, const vector<glm::vec4> &segments) {
	Scalar rnd_color = generate_random_color();
	draw_segments_on_image(img, segments,rnd_color);
}

void draw_segments_on_image_rnd_colors(Mat &img, const vector<glm::vec4> &segments) {
	vector<Scalar> rnd_colors = generate_random_colors(segments.size());
	draw_segments_on_image(img, segments,rnd_colors);
}

void draw_segments_on_image(Mat &img, const vector<glm::vec4> &segments, const Scalar color) {
	for(const auto &s : segments)
		draw_segment_on_img(img,s,color);
}

void draw_segments_on_image(Mat &img, const vector<glm::vec4> &segments, const Scalar color, float thickness) {
	for(const auto &s : segments)
		draw_segment_on_img(img,s,color,thickness);
}

void draw_segments_on_image(Mat &img, const vector<glm::vec4> &segments, const vector<Scalar> colors) {
	for(ulong i=0; i < segments.size();i++)
		draw_segment_on_img(img,segments[i],colors[i]);
}

void draw_segments_on_image(Mat &img, const vector<Vec4f> &segments) {
	Ptr<LineSegmentDetector> ls = createLineSegmentDetector(LSD_REFINE_ADV);
	ls->drawSegments(img, segments);
}

Mat get_black_image(const Mat &img) {
	Mat bim = Mat(img.rows,img.cols,CV_8UC3);
	bim.setTo(cv::Scalar(0,0,0));
	return bim;
}

/**
 * Build new square image centered around coords from given img, with given 2*range+1
 */
Mat getFocusImage(const Mat &img, const glm::vec2 &coords, const float range, const int image_type) {
	int square_dim = range*2+1;

	const int miny = max(0,int(ceil(coords.y-range)));
	const int maxy = min(img.rows,miny+square_dim-1);

	const int minx = max(0,int(ceil(coords.x-range)));
	const int maxx = min(img.cols,minx+square_dim-1);

	Mat new_img = Mat(square_dim,square_dim,image_type);
	img(Range(miny,maxy),Range(minx,maxx)).copyTo(new_img);

	return new_img;
}

/**
 * Given a vector of imgs, choose focus points using focus_points_img_ids and focus_points_coords.
 * Build new square images centered around focus points, with given 2*range+1, and return them
 * as a vector
 */
vector<Mat> getVectorFocusImages(const vector<Mat> &imgs,const vector<int> &focus_points_img_ids, const vector<glm::vec2> &focus_points_coords, const int range, const int image_type) {
	vector<Mat> new_imgs;
	int img_id;
	glm::vec2 coords;


	for(int i=0;i<focus_points_img_ids.size();i++) {
		img_id = focus_points_img_ids[i];
		const Mat &img = imgs[img_id];
		coords = focus_points_coords[i];
		new_imgs.push_back(getFocusImage(img, coords, range,image_type));
	}

	return new_imgs;
}

/**
 * Given a vector of imgs, choose focus points using focus_points_img_ids and focus_points_coords.
 * Build new square images centered around focus points, with given 2*range+1, and return them
 * as a unified focus image
 */
Mat getUnifiedFocusImages(const vector<Mat> &imgs,const vector<int> &focus_points_img_ids, const vector<glm::vec2> &focus_points_coords, const int range, const int image_type) {
	vector<Mat> focus_images = getVectorFocusImages(imgs,focus_points_img_ids, focus_points_coords, range,image_type);

	Mat focus_image = getUnifiedSquareImage(focus_images, image_type);

	release_vector_Mat(focus_images);

	return focus_image;
}

void draw_single_point_process_with_segment_output_nosegmdrawing(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<SegmentEdgeManager::segment_point> &starting_intersections, const int correspondence_detection_radius, const std::vector< std::vector<std::vector<SegmentEdgeManager::segment_point>>> &detected_correspondent_intersections) {
	Scalar starting_refpoint_color = generate_random_color();
	vector<Scalar> starting_intersection_colors = generate_random_colors(starting_intersections.size());
	draw_single_point_process_with_segment_output_nosegmdrawing(imgs, sfmd, all_fundamental_matrices, starting_img, reference_point_id, initial_detection_radius, starting_intersections, correspondence_detection_radius, detected_correspondent_intersections, starting_refpoint_color ,starting_intersection_colors);
}

/*
* detected_correspondent_intersections : for each starting intersection point, for each image where reference point is visible except starting image, each correspondent intersection found
*/
void draw_single_point_process_with_segment_output_nosegmdrawing(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<SegmentEdgeManager::segment_point> &starting_intersections, const int correspondence_detection_radius, const std::vector< std::vector<std::vector<SegmentEdgeManager::segment_point>>> &detected_correspondent_intersections, const Scalar &starting_refpoint_color ,const vector<Scalar> &starting_intersection_colors) {
	glm::vec2 p_cur2dcoords;
	const std::vector<int> &cams = sfmd.camViewingPointN_[reference_point_id];
	const std::vector<glm::vec2> &coords_on_cams = sfmd.point2DoncamViewingPoint_[reference_point_id];
	int cur_cam;
	int cam_index,starting_intersection_index;
	const int starting_intersections_num = starting_intersections.size();
	glm::vec3 epipolar_line;

	glm::vec2 refpoint_starting_image;
	bool found;
	get_2d_coordinates_of_point_on_image(sfmd,starting_img,reference_point_id,refpoint_starting_image,found);


	// Draw reference point and, only on starting image , corresponding detection radius and starting intersections
	for(cam_index=0; cam_index < cams.size(); cam_index++) {
		cur_cam = cams[cam_index];
		p_cur2dcoords = coords_on_cams[cam_index];

		// Draw reference starting point on each image where it's visible
		draw_reference_point_glm(imgs[cur_cam], p_cur2dcoords, starting_refpoint_color);

		if(cur_cam == starting_img) {
			// Draw starting detection radius
			draw_circle_glm(imgs[starting_img], p_cur2dcoords,initial_detection_radius, starting_refpoint_color);

			for(starting_intersection_index=0; starting_intersection_index < starting_intersections_num; starting_intersection_index++) {
				// Draw starting intersection on first image
				draw_intersection_point_glm(imgs[starting_img], starting_intersections[starting_intersection_index].coords, starting_intersection_colors[starting_intersection_index]);
			}
		} else {
			// Draw max range detection radius on other images
			draw_circle_glm(imgs[cur_cam], p_cur2dcoords,correspondence_detection_radius, starting_refpoint_color);

			// Get Fundamental Matrix
			cv::Mat &F = all_fundamental_matrices[starting_img][cur_cam];

			// Draw epipolar line on refpoint
			if(computeCorrespondEpilineSinglePoint(refpoint_starting_image, F, epipolar_line,1))
				draw_line_glm(imgs[cur_cam], epipolar_line, starting_refpoint_color);

			for(starting_intersection_index=0; starting_intersection_index < starting_intersections_num; starting_intersection_index++) {
				//Draw epipolar line
				if(computeCorrespondEpilineSinglePoint(starting_intersections[starting_intersection_index].coords, F, epipolar_line,1))
					draw_line_glm(imgs[cur_cam], epipolar_line, starting_intersection_colors[starting_intersection_index]);

				for(const auto &correspondence_point : detected_correspondent_intersections[starting_intersection_index][cam_index]) {
					// Draw starting intersection on first image
					draw_intersection_point_glm(imgs[cur_cam], correspondence_point.coords, starting_intersection_colors[starting_intersection_index]);
				}
			}
		}
	}
}

Mat convert_to_lsd_edgesimgs(const Mat &img) {
	vector<Vec4f> cur_edges;
	Ptr<LineSegmentDetector> ls = createLineSegmentDetector(LSD_REFINE_ADV);
	Mat cur_gs;
	cvtColor( img, cur_gs, CV_BGR2GRAY );
	ls->detect(cur_gs, cur_edges);
	Mat ni = get_black_image(img);
	draw_segments_on_image(ni,cur_edges);
	return ni;
}

vector<Mat> convert_to_lsd_edgesimgs(const vector<Mat> &imgs) {
	vector<Mat> nis;

	for(const auto &img : imgs)
		nis.push_back(convert_to_lsd_edgesimgs(img));

	return nis;
}

Mat draw_MultiColorComponents_PolyLineGraph_simplified(const Mat &img, const PolyLineGraph2D & plg) {
	vector<vector<glm::vec4>> vec_segments = plg.get_segments_grouped_by_component();
	Mat outimg = get_black_image(img);
	for(const auto &segments : vec_segments)
		draw_segments_on_image_rnd_color(outimg,segments);

	draw_points_glm(outimg, plg.get_nodes_with_loops_coords(), Scalar(255,255,255), 2);

	draw_points_glm(outimg, plg.get_hub_nodes_coords(), Scalar(255,0,0), 1);
	draw_points_glm(outimg, plg.get_extreme_nodes_coords(), Scalar(0,255,0), 1);
	draw_points_glm(outimg, plg.get_loopnodes_coords(), Scalar(0,0,255), 1);

	return outimg;
}

Mat draw_MultiColorSegments_PolyLineGraph_simplified(const Mat &img, const PolyLineGraph2D & plg) {
	vector<glm::vec4> segments = plg.get_segments_list();
	Mat outimg = get_black_image(img);
	draw_segments_on_image_rnd_colors(outimg,segments);

	vector<glm::vec2> segments_extremes;
	for(const auto &s : segments) {
		segments_extremes.push_back(glm::vec2(s[0],s[1]));
		segments_extremes.push_back(glm::vec2(s[2],s[3]));
	}

	draw_points_glm(outimg, segments_extremes, Scalar(255,255,255), 1);

	return outimg;
}

vector<Mat> draw_plgs(const vector<Mat> &imgs, const vector<PolyLineGraph2DHMapImpl> &plgs) {
	vector<Mat> res;

	for(int i=0; i < imgs.size(); i++) {
		const PolyLineGraph2DHMapImpl &plg = plgs[i];
		const Mat &img = imgs[i];
		res.push_back(draw_MultiColorComponents_PolyLineGraph_simplified(img, plg));
	}

	return res;
}

vector<Mat> draw_MultiColorSegments_PolyLineGraphs_simplified(const vector<Mat> &imgs, const vector<PolyLineGraph2DHMapImpl> &plgs) {
	vector<Mat> res;

	for(int i=0; i < imgs.size(); i++) {
		const PolyLineGraph2DHMapImpl &plg = plgs[i];
		const Mat &img = imgs[i];
		res.push_back(draw_MultiColorSegments_PolyLineGraph_simplified(img, plg));
	}

	return res;
}

vector<Mat> draw_MultiColorPolyLines_PolyLineGraphs_simplified(const vector<Mat> &imgs, const vector<PolyLineGraph2DHMapImpl> &plgs) {
	vector<Mat> res;

	for(int i=0; i < imgs.size(); i++) {
		const PolyLineGraph2DHMapImpl &plg = plgs[i];
		const Mat &img = imgs[i];
		res.push_back(draw_MultiColorPolyLines_PolyLineGraph_simplified(img, plg));
	}

	return res;
}

Mat draw_overlay_MultiColorComponents_PolyLineGraph_simplified(Mat &img, const PolyLineGraph2D & plg) {
	vector<vector<glm::vec4>> vec_segments = plg.get_segments_grouped_by_component();
	//vector<glm::vec4> segments = convertEdgeImagePixelToSegments_PolyLineGraph(img,edge_color);
	cout << "PLG Connected components: " << vec_segments.size() << endl;
	for(const auto &segments : vec_segments)
		draw_segments_on_image_rnd_color(img,segments);

	draw_points_glm(img, plg.get_nodes_with_loops_coords(), Scalar(255,255,255), 2);

	draw_points_glm(img, plg.get_hub_nodes_coords(), Scalar(255,0,0), 1);
	draw_points_glm(img, plg.get_extreme_nodes_coords(), Scalar(0,255,0), 1);
	draw_points_glm(img, plg.get_loopnodes_coords(), Scalar(0,0,255), 1);

	return img;
}

Mat draw_MultiColorPolyLines_PolyLineGraph_simplified(const Mat &img, const PolyLineGraph2D & plg) {
	vector<vector<glm::vec4>> vec_segments = plg.get_segments_grouped_by_polyline();
	Mat outimg = get_black_image(img);
	for(const auto &segments : vec_segments)
		draw_segments_on_image_rnd_color(outimg,segments);
	return outimg;
}

Mat draw_polyline_graph_simplified(const Mat &img, const PolyLineGraph2D & plg, const Scalar &color) {
	Mat outimg = get_black_image(img);
	draw_PolyLineGraph_simplified_overlay(outimg, plg, color); 
	return outimg;
}

void draw_PolyLineGraph_simplified_overlay(Mat &img, const PolyLineGraph2D & plg, const Scalar &color) {
	vector<glm::vec4> segments = plg.get_segments_list();
	draw_segments_on_image(img,segments,color);
}

Mat create_image_frame(const Mat &img, const Scalar &frameColor, const int frameWidth) {
	Mat out = Mat(img.rows+2*frameWidth,img.cols+2*frameWidth,CV_8UC3);

	for(int y=0; y < img.rows; y++)
		for(int x=0; x < img.cols; x++)
			out.at<Vec3b>(y+frameWidth,x+frameWidth) = img.at<Vec3b>(y,x);

	for(int y=0; y < frameWidth; y++)
		for(int x=0; x < img.cols+2*frameWidth; x++)
			out.at<Vec3b>(y,x) = Vec3b(0,0,0);

	for(int y=0; y < frameWidth; y++)
		for(int x=0; x < img.cols+2*frameWidth; x++)
			out.at<Vec3b>(out.rows-1-y,x) = Vec3b(0,0,0);

	for(int y=0; y < img.rows+2*frameWidth; y++)
		for(int x=0; x < frameWidth; x++)
			out.at<Vec3b>(y,x) = Vec3b(0,0,0);

	for(int y=0; y < img.rows+2*frameWidth; y++)
		for(int x=0; x < frameWidth; x++)
			out.at<Vec3b>(y,out.cols-1-x) = Vec3b(0,0,0);

	return out;
}

void draw_colored_polylines_lsd(const vector<Mat> &imgs, const SfMData &sfmd ,char *em_out_folder) {
	vector<Mat> plgs_imgs =  convertEdgeImagesPixelToSegmentsImages_PolyLineGraph_simplified(imgs,EDGE_COLOR);
	vector<Mat> lsd_plgs_imgs = convert_to_lsd_edgesimgs(plgs_imgs);
	write_images(em_out_folder, sfmd, lsd_plgs_imgs, "polylines_lsd_");
}

void draw_colored_components_and_edge_refpoints(const vector<Mat> &imgs, const SfMData &sfmd ,char *em_out_folder) {
	vector<PolyLineGraph2DHMapImpl> plgs = convertEdgeImagesPolyLineGraphs_simplified(imgs,EDGE_COLOR);
	std::vector<int> edgerefpoints = find_edgerefpoints(plgs,sfmd);
	vector<Mat> out_imgs;
	for(int i=0; i < plgs.size(); i++) {
		const PolyLineGraph2DHMapImpl &plg = plgs[i];
		const Mat &img = imgs[i];
		out_imgs.push_back(draw_MultiColorComponents_PolyLineGraph_simplified(img,plg));
	}
	draw_setofrefpoints_on_imgs(out_imgs,sfmd, edgerefpoints);
	write_images(em_out_folder, sfmd, out_imgs, "segmented_input_edges_");
}

void draw_polyline_matches(vector<Mat> &imgs, const vector<PolyLineGraph2DHMapImpl> &plgs, const vector<vector<set<ulong>>> &pl_matches)
{
	for(const auto &pl_match : pl_matches) {
		const Scalar col = generate_random_color();

		for(int i=0; i < plgs.size(); i++)
			for(set<ulong>::iterator jt=pl_match[i].begin(); jt!=pl_match[i].end(); jt++)
				draw_segments_on_image(imgs[i], plgs[i].polylines[*jt].get_segments_list(),col);
	}
}

void draw_and_write_focus_image(vector<Mat> &imgs, const SfMData &sfm_data, const int refpoint_id, const int img_to_start_id, int global_processed_point_counter, char * outfolder) {
	Mat unifiedFocusImage = getUnifiedFocusImages(imgs,sfm_data.camViewingPointN_[refpoint_id],sfm_data.point2DoncamViewingPoint_[refpoint_id], FOCUS_IMAGES_RANGE, CV_8UC3);
	string img_path = outfolder + std::string("_focuspoint") + std::to_string(global_processed_point_counter) + std::string("_starting") + std::to_string(img_to_start_id) + std::string("_refpoint") + std::to_string(refpoint_id) +".png";
	imwrite(img_path,unifiedFocusImage);
	unifiedFocusImage.release();
	global_processed_point_counter++;
}

void draw_sfmd_points_overwrite(vector<Mat> &out_imgs, const vector<PolyLineGraph2DHMapImpl> &plgs, const SfMData &sfm_data, const ulong start_from_point, const char *out_folder, const String &s) {
	for(int point_id=start_from_point; point_id < sfm_data.numPoints_; point_id++)
		draw_refpoint_on_imgs(out_imgs,sfm_data, point_id, generate_random_color());

	write_images(out_folder, sfm_data, out_imgs, s);
}

vector<Mat> draw_plgs_bw(vector<Mat> &imgs, const vector<PolyLineGraph2DHMapImpl> &plgs) {
	vector<Mat> out_imgs;
	for(int i=0; i < plgs.size(); i++)
		out_imgs.push_back(draw_polyline_graph_simplified(imgs[i],plgs[i],Scalar(255,255,255)));
	return out_imgs;
}

void draw_sfmd_points(vector<Mat> &imgs, const vector<PolyLineGraph2DHMapImpl> &plgs, const SfMData &sfm_data, const char *out_folder, const String &s) {
	draw_sfmd_points(imgs, plgs, sfm_data, 0, out_folder, s);
}

void draw_sfmd_points(vector<Mat> &imgs, const vector<PolyLineGraph2DHMapImpl> &plgs, const SfMData &sfm_data, const ulong start_from_point, const char *out_folder, const String &s) {
	vector<Mat> out_imgs_orig_bg = copy_vector_Mat(imgs);
	for(int i=0; i < plgs.size(); i++)
		draw_PolyLineGraph_simplified_overlay(out_imgs_orig_bg[i],plgs[i],Scalar(0,0,0));
	draw_sfmd_points_overwrite(out_imgs_orig_bg, plgs, sfm_data, start_from_point, out_folder, s);
}

void draw_sfmd_points_plgs(vector<Mat> &imgs, const vector<PolyLineGraph2DHMapImpl> &plgs, const SfMData &sfm_data, const char *out_folder, const String &s) {
	draw_sfmd_points_plgs(imgs, plgs, sfm_data, 0, out_folder, s);
}

void draw_sfmd_points_plgs(vector<Mat> &imgs, const vector<PolyLineGraph2DHMapImpl> &plgs, const SfMData &sfm_data, const ulong start_from_point, const char *out_folder, const String &s) {
	vector<Mat> out_imgs_plgs = draw_plgs_bw(imgs, plgs);
	draw_sfmd_points_overwrite(out_imgs_plgs, plgs, sfm_data, start_from_point, out_folder, s);
}




