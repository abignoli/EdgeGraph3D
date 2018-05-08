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


#ifndef INCLUDE_EDGEGRAPH3D_UTILS_DRAWING_UTILITIES_HPP_
#define INCLUDE_EDGEGRAPH3D_UTILS_DRAWING_UTILITIES_HPP_

#include <opencv2/core/mat.hpp>
#include <opencv2/core/matx.hpp>
#include <opencv2/core/types.hpp>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "segment_edge_manager.hpp"
#include "polyline_graph_2d.hpp"
#include "polyline_graph_2d_hmap_impl.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"
#include "glm.hpp"


using namespace std;
using namespace cv;

using namespace cv;

#define DRAW_EMPTY_CIRCLE 1
#define DRAW_FILLED_CIRCLE -1
#define DRAW_REFERENCE_POINT_RADIUS 2
#define DRAW_INTERSECTION_POINT_RADIUS 2
#define DRAW_NEW_MATCHED_POINT_RADIUS (DRAW_INTERSECTION_POINT_RADIUS+1)
#define DRAW_CORRESPONDENCE_INNER (DRAW_INTERSECTION_POINT_RADIUS-1 > 2 ? DRAW_INTERSECTION_POINT_RADIUS-1 : 2)

void draw_point(Mat &img, const Point2f &p, const Scalar &color, const int draw_radius);

void draw_point_glm(Mat &img, const glm::vec2 &p, const Scalar &color, const int draw_radius);

void draw_point(Mat &img, const Point2f &p, const int draw_radius);

void draw_point_glm(Mat &img, const glm::vec2 &p, const int draw_radius);

void draw_points_glm(Mat &img, const vector<glm::vec2> &ps, const int draw_radius);

void draw_points_glm(Mat &img, const vector<glm::vec2> &ps, const Scalar &color, const int draw_radius);

void draw_reference_point_glm(Mat &img, const glm::vec2 &p, const Scalar &color);

void draw_intersection_point_glm(Mat &img, const glm::vec2 &p, const Scalar &color);

void draw_line_glm(Mat& img, const glm::vec3 &line, const Scalar& color);

void draw_circle_glm(Mat &img, const glm::vec2 &center,const int radius, const Scalar &color);

void draw_single_point_process(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<glm::vec2> &starting_intersections, const int correspondence_detection_radius, const std::vector< std::vector<std::vector<glm::vec2>>> &detected_correspondent_intersections, const Scalar &starting_refpoint_color ,const vector<Scalar> &starting_intersection_colors);

void draw_single_point_process(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<glm::vec2> &starting_intersections, const int correspondence_detection_radius, const std::vector< std::vector<std::vector<glm::vec2>>> &detected_correspondent_intersections);

void draw_single_point_process(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<PolyLineGraph2D::plg_point> &starting_intersections_plgp, const std::vector< std::vector<std::vector<PolyLineGraph2D::plg_point>>> &detected_correspondent_intersections_plgp);

void draw_single_point_process(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<pair<PolyLineGraph2D::plg_point,std::vector<std::vector<PolyLineGraph2D::plg_point>>>> &nearby_matches);

void draw_single_point_process_no_epilines(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<pair<PolyLineGraph2D::plg_point,std::vector<std::vector<PolyLineGraph2D::plg_point>>>> &nearby_matches);

void draw_refpoint_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const int point_id, const Scalar &point_color);

void draw_refpoints_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const vector<Scalar> &point_colors);

void draw_refpoint_with_circle_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const int point_id, const Scalar &point_color, const float radius);

void draw_refpoints_with_circle_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const vector<Scalar> &point_colors, const float radius);

void draw_refpoint_with_circles_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const int point_id, const Scalar &point_color, const float radius1, const float radius2);

void draw_refpoints_with_circles_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const vector<Scalar> &point_colors, const float radius1, const float radius2);

void draw_setofrefpoints_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const vector<int> refpoints, const vector<Scalar> &point_colors);

void draw_setofrefpoints_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const vector<int> refpoints);

void draw_refpoint_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const int point_id);

void draw_refpoints_on_imgs(vector<Mat> &imgs,const SfMData &sfmd);

void draw_point_projections(vector<Mat> &imgs, const vector<glm::vec2> &coords, const vector<int> &cameras, const Scalar &point_color);
void draw_point_projections(vector<Mat> &imgs, const vector<glm::vec2> &coords, const vector<int> &cameras);

void draw_3dpoint_on_imgs(vector<Mat> &imgs, const std::tuple<glm::vec3, vector<glm::vec2>, vector<int>> &p3d, const Scalar &point_color);
void draw_3dpoint_on_imgs(vector<Mat> &imgs, const std::tuple<glm::vec3, vector<glm::vec2>, vector<int>> &p3d);

void draw_3dpoints_on_imgs(vector<Mat> &imgs, const vector<std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>> &p3ds);

void draw_point_projections(vector<Mat> &imgs, const vector<PolyLineGraph2D::plg_point> &coords, const vector<int> &cameras, const Scalar &point_color);

void draw_point_projections(vector<Mat> &imgs, const vector<PolyLineGraph2D::plg_point> &coords, const vector<int> &cameras);

void draw_3dpoint_on_imgs(vector<Mat> &imgs, const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &p3d, const Scalar &point_color);

void draw_3dpoint_on_imgs(vector<Mat> &imgs, const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &p3d);

void draw_3dpoints_on_imgs(vector<Mat> &imgs, const vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &p3ds);

void draw_refpoint_with_circle_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const int point_id, const float radius);

void draw_refpoints_with_circle_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const float radius);

void draw_refpoint_with_circles_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const int point_id, const float radius1, const float radius2);

void draw_refpoints_with_circles_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, const float radius1, const float radius2);

void draw_point_epipolars_on_imgs(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, glm::vec2 p, const int starting_img, const Scalar &point_color);

void draw_refpoint_epipolars_on_imgs(vector<Mat> &imgs,const SfMData &sfmd,const cv::Mat** all_fundamental_matrices, const int point_id, const int starting_img, const Scalar &point_color);

void draw_img_pair_refpoints(vector<Mat> &imgs,const SfMData &sfmd,const int i,const int j);

void draw_img_pair_epipolars_refpoints(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices,const int i,const int j);

char showImages(string title,const vector<Mat>& imgs,const Size &totalSize, const int image_type);

Scalar generate_random_color();

vector<Scalar> generate_random_colors(int size);

void draw_consensus_matched_point_on_img(const Mat &img,const glm::vec2 &coords, const Scalar &color);

void draw_consensus_matched_point(const vector<Mat>& imgs,const vector<int> &cam_ids, const vector<glm::vec2> &cam_coords, const Scalar &color);

void draw_consensus_matched_point(const vector<Mat>& imgs,const vector<int> &cam_ids, const vector<glm::vec2> &cam_coords);

void draw_consensus_matched_points(const vector<Mat>& imgs,const vector<std::tuple<glm::vec3,vector<glm::vec2>,vector<int>>> &new_points);

void draw_consensus_matched_points(const vector<Mat>& imgs,const vector<std::tuple<glm::vec3,vector<glm::vec2>,vector<int>>> &new_points, const Scalar &color);

void draw_single_point_process_and_consensus(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<glm::vec2> &starting_intersections, const int correspondence_detection_radius, const std::vector< std::vector<std::vector<glm::vec2>>> &detected_correspondent_intersections, const vector<std::tuple<glm::vec3,vector<glm::vec2>,vector<int>>> &found_3d_points, const Scalar &consensus_match_color);

void draw_single_point_process_and_consensus(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<SegmentEdgeManager::segment_point> &starting_intersections_w_segments, const int correspondence_detection_radius, const std::vector< std::vector<std::vector<SegmentEdgeManager::segment_point>>> &detected_correspondent_intersections_w_segments, const vector<std::tuple<glm::vec3,vector<glm::vec2>,vector<int>>> &found_3d_points, const Scalar &consensus_match_color);

void draw_new_consensus_points(vector<Mat> &imgs, const vector<std::tuple<glm::vec3,vector<glm::vec2>,vector<int>>> &found_3d_points);

void draw_segments_on_image(Mat &img, const vector<glm::vec4> &segments);

//void draw_segments_on_image(Mat &out, const vector<glm::vec4> &segments, const Scalar &colorlines);

Mat draw_segments_on_newimage(const Size &sz, const vector<glm::vec4> &segments, const Scalar &colorbg, const Scalar &colorlines);

Mat draw_segments_on_newimage_with_extremes(const Size &sz, const vector<glm::vec4> &segments, const Scalar &colorbg, const Scalar &colorlines, const Scalar &colorstart, const Scalar &colorend);

void draw_segments_on_image(Mat &img, const vector<Vec4f> &segments);

Mat get_black_image(const Mat &img);

Mat draw_MultiColorComponents_PolyLineGraph_simplified(const Mat &img, const PolyLineGraph2D & plg);

Mat draw_MultiColorSegments_PolyLineGraph_simplified(const Mat &img, const PolyLineGraph2D & plg);

vector<Mat> draw_plgs(const vector<Mat> &imgs, const vector<PolyLineGraph2DHMapImpl> &plgs);

vector<Mat> draw_MultiColorSegments_PolyLineGraphs_simplified(const vector<Mat> &imgs, const vector<PolyLineGraph2DHMapImpl> &plgs);

vector<Mat> draw_MultiColorPolyLines_PolyLineGraphs_simplified(const vector<Mat> &imgs, const vector<PolyLineGraph2DHMapImpl> &plgs);

Mat draw_overlay_MultiColorComponents_PolyLineGraph_simplified(Mat &img, const PolyLineGraph2D & plg);

/**
 * Given a vector of imgs, choose focus points using focus_points_img_ids and focus_points_coords.
 * Build new square images centered around focus points, with given 2*range+1, and return them
 * as a unified focus image
 */
Mat getUnifiedFocusImages(const vector<Mat> &imgs,const vector<int> &focus_points_img_ids, const vector<glm::vec2> &focus_points_coords, const int range, const int image_type);

/**
 * Images are arranged in a squared display, with no scaling.
 */
Mat getUnifiedSquareImage(const vector<Mat>& imgs, const int image_type);

/**
 * Images are arranged in a squared display, of the specified size.
 */
Mat getUnifiedSquareImage(const vector<Mat>& imgs,const Size &totalSize, const int image_type);

void draw_single_point_process_with_segment_output(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<SegmentEdgeManager::segment_point> &starting_intersections, const int correspondence_detection_radius, const std::vector< std::vector<std::vector<SegmentEdgeManager::segment_point>>> &detected_correspondent_intersections);

void draw_single_point_process_with_segment_output(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<SegmentEdgeManager::segment_point> &starting_intersections, const int correspondence_detection_radius, const std::vector< std::vector<std::vector<SegmentEdgeManager::segment_point>>> &detected_correspondent_intersections, const Scalar &starting_refpoint_color ,const vector<Scalar> &starting_intersection_colors);

void draw_single_point_process_with_segment_output_nosegmdrawing(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<SegmentEdgeManager::segment_point> &starting_intersections, const int correspondence_detection_radius, const std::vector< std::vector<std::vector<SegmentEdgeManager::segment_point>>> &detected_correspondent_intersections);

void draw_single_point_process_with_segment_output_nosegmdrawing(vector<Mat> &imgs,const SfMData &sfmd, cv::Mat** all_fundamental_matrices, const int starting_img, const int reference_point_id, const int initial_detection_radius, const std::vector<SegmentEdgeManager::segment_point> &starting_intersections, const int correspondence_detection_radius, const std::vector< std::vector<std::vector<SegmentEdgeManager::segment_point>>> &detected_correspondent_intersections, const Scalar &starting_refpoint_color ,const vector<Scalar> &starting_intersection_colors);

void draw_segments_on_image_rnd_color(Mat &img, const vector<glm::vec4> &segments);

void draw_segments_on_image_rnd_colors(Mat &img, const vector<glm::vec4> &segments);

void draw_segments_on_image(Mat &img, const vector<glm::vec4> &segments, const Scalar color);
void draw_segments_on_image(Mat &img, const vector<glm::vec4> &segments, const Scalar color, float thickness);

void draw_segments_on_image(Mat &img, const vector<glm::vec4> &segments, const vector<Scalar> colors);

Mat convert_to_lsd_edgesimgs(const Mat &img);

vector<Mat> convert_to_lsd_edgesimgs(const vector<Mat> &imgs);

Mat draw_MultiColorPolyLines_PolyLineGraph_simplified(const Mat &img, const PolyLineGraph2D & plg);

Mat draw_polyline_graph_simplified(const Mat &img, const PolyLineGraph2D & plg, const Scalar &color);

void draw_PolyLineGraph_simplified_overlay(Mat &img, const PolyLineGraph2D & plg, const Scalar &color);

Mat create_image_frame(const Mat &img, const Scalar &frameColor, const int frameWidth);

void draw_colored_polylines_lsd(const vector<Mat> &imgs, const SfMData &sfmd ,char *em_out_folder);

void draw_colored_components_and_edge_refpoints(const vector<Mat> &imgs, const SfMData &sfmd ,char *em_out_folder);

void draw_polyline_matches(vector<Mat> &imgs, const vector<PolyLineGraph2DHMapImpl> &plgs, const vector<vector<set<ulong>>> &pl_matches);

void draw_and_write_focus_image(vector<Mat> &imgs, const SfMData &sfm_data, const int refpoint_id, const int img_to_start_id, int global_processed_point_counter, char * outfolder);

void draw_sfmd_points_overwrite(vector<Mat> &out_imgs, const vector<PolyLineGraph2DHMapImpl> &plgs, const SfMData &sfm_data, const ulong start_from_point, const char *out_folder, const String &s);

vector<Mat> draw_plgs_bw(vector<Mat> &imgs, const vector<PolyLineGraph2DHMapImpl> &plgs);

void draw_sfmd_points(vector<Mat> &imgs, const vector<PolyLineGraph2DHMapImpl> &plgs, const char *out_folder, const String &s);

void draw_sfmd_points(vector<Mat> &imgs, const vector<PolyLineGraph2DHMapImpl> &plgs, const SfMData &sfm_data, const ulong start_from_point, const char *out_folder, const String &s);

void draw_sfmd_points_plgs(vector<Mat> &imgs, const vector<PolyLineGraph2DHMapImpl> &plgs, const SfMData &sfm_data, const char *out_folder, const String &s);

void draw_sfmd_points_plgs(vector<Mat> &imgs, const vector<PolyLineGraph2DHMapImpl> &plgs, const SfMData &sfm_data, const ulong start_from_point, const char *out_folder, const String &s);

#endif /* INCLUDE_EDGEGRAPH3D_UTILS_DRAWING_UTILITIES_HPP_ */
