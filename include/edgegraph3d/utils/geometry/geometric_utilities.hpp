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


#ifndef INCLUDE_EDGEGRAPH3D_UTILS_GEOMETRY_GEOMETRIC_UTILITIES_HPP_
#define INCLUDE_EDGEGRAPH3D_UTILS_GEOMETRY_GEOMETRIC_UTILITIES_HPP_

#include <CGAL/Cartesian.h>
#include <CGAL/Line_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>

#include <utility>
#include <vector>

#include "segment_edge_manager.hpp"
#include "datatypes.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"
#include "glm.hpp"


struct CameraType;

#define MAX_LINE3DANGLE_COS 0.5

using namespace std;
using namespace cv;

typedef CGAL::Point_3<CGAL::Cartesian<double>> CGAL_Point3;
typedef CGAL::Line_3<CGAL::Cartesian<double>> CGAL_Line3;
typedef CGAL::Vector_3<CGAL::Cartesian<double>> CGAL_Vector3;
typedef CGAL::Plane_3<CGAL::Cartesian<double>> CGAL_Plane3;

/**
 * Output: [a,b,c,d]
 * Plane is in the form: a * x + b * y + c * z + d = 0
 */
glm::vec4 plane_from_3_points_m(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3);

std::pair<glm::vec3,float> plane_from_3_points_pairoutput_m(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3);

/**
 * Output: p1, v1, v2
 * Plane is in the form: p1 + k * v1 + t * v2
 */
glm::mat3 plane_from_3_points_2vecs_form_m(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3);

/**
 * Computes the intersection between two planes:
 * 	- if parallel, set parallel true
 * 	- if overlapped (i.e. same plane), set overlapped true
 * 	- otherwise define the intersection line as:
 * 		l : x0 + k * v, k in R
 *		returned as line3d = (x0,v)
 *
 */
void planes_intersection_m(const std::pair<glm::vec3,float> &p1, const std::pair<glm::vec3,float> &p2, std::pair<glm::vec3,glm::vec3> &line3d, bool &parallel, bool &overlapped);


/**
 * Compute Fundamental matrices for each image pair
 */
cv::Mat** generate_all_fundamental_matrices(const SfMData &sfmd);


/**
 * Computes intersections between line described by (a,b,c) and the circle of radius r centered in (xc,yc).
 *
 * Returns vector of points of length:
 * 	0 -> no intersections
 * 	1 -> 1 intersection (i.e. line is tangent)
 * 	2 -> 2 intersections
 */
vector<glm::vec2> detect_circle_line_intersections(const glm::vec3 &line_params,const glm::vec2 &center,const float r);

std::pair<CGAL_Plane3, std::vector<CGAL_Plane3>> compute_3d_planes_potential_match(const SfMData &sfmd, const int starting_img_id, const int refpoint, const potential_match &pm);

/**
 * Computes intersections between segment from (x1,y1) to (x2,y2) and circle of radius r centered in (xc,yc)
 */
vector<glm::vec2> detect_circle_segment_intersections(const glm::vec4 &segm,const glm::vec2 &center,const float r);

/*template <typename T>
void FreeAll( T & t ) {
    T tmp;
    t.swap( tmp );
}*/

/**
 * Compute intersection between segment and line. Possible results are:
 * 	- parallel
 * 		- overlapped
 * 		- not overlapped
 * 	- single intersection
 *
 * 	Usage should first check parallel and overlapped flags, before reading the intersection field.
 */
void intersect_segment_line(const glm::vec4 &segm,const glm::vec3 &line, bool &parallel, bool &overlapped, bool &intersection_found, glm::vec2 &intersection);

glm::vec3 get_line_from_ray(const ray2d &r);

void intersect_segment_ray(const glm::vec4 &segm,const ray2d &r, bool &parallel, bool &overlapped, bool &past_center, bool &intersection_found, glm::vec2 &intersection);

bool aligned(const glm::vec2 &a,const glm::vec2 &b,const glm::vec2 &c);

/**
 * Compute intersection between segment and line. Possible results are:
 * 	- parallel
 * 		- overlapped
 * 		- not overlapped
 * 	- single intersection
 * 	- quasiparallel
 * 		- if direction is almost parallel: cos > max_quasiparallel_angle_cos
 * 		- if distance of closest point < max_quasiparallel_dist
 * 	- valid
 * 		- true : if no intersection is found, or it is found and the segment is not quasiparallel within max_dist
 * 		- false : segment is quasiparallel within max_dist (intersection can be found or not)
 *
 * 	valid flag indicates where the intersection computation is considered accurate enough. If line and segment
 * 	are quasiparallel, the intersection point (or the non-intersection) cannot be established in a sufficiently
 * 	precise manner
 *
 * 	Usage should first check parallel and overlapped flags, before reading the intersection field.
 */
void intersect_segment_line_no_quasiparallel(const glm::vec4 &segm,const glm::vec3 &line, const float max_quasiparallel_angle_cos, const float max_quasiparallel_dist, bool &parallel, bool &overlapped, bool &intersection_found, bool &quasiparallel_within_distance, bool &valid, float &distance, glm::vec2 &intersection);
bool point_in_segment_bounding_box(const glm::vec4 &segm, const glm::vec2 &pt);

void intersect_segment_segment(const glm::vec4 &segm1,const glm::vec4 &segm2, bool &parallel, bool &overlapped, bool &intersection_found, glm::vec2 &intersection);

void intersect_line_line(const glm::vec3 &line_i,const glm::vec3 &line_j, bool &parallel, bool &overlapped, glm::vec2 &intersection);

float squared_2d_distance(const glm::vec2 &a,const glm::vec2 &b);

float squared_3d_distance(const glm::vec3 &a,const glm::vec3 &b);

float compute_lengthsq(const glm::vec2 &a);
float compute_lengthsq(const glm::vec3 &a);

float compute_2d_distance(const glm::vec2 &a,const glm::vec2 &b);

float compute_3d_distance(const glm::vec3 &a,const glm::vec3 &b);

float compute_anglecos_vec2_vec2(const glm::vec2 &a,const glm::vec2 &b);

float compute_anglecos_vec2_vec2(const glm::vec2 &a1,const glm::vec2 &a2,const glm::vec2 &b1,const glm::vec2 &b2);

glm::vec2 get_2d_direction(const glm::vec4 &segm);

glm::vec2 get_2d_direction(const glm::vec3 &line);

float compute_anglecos(const glm::vec4 &segm,const glm::vec3 &line);

/**
 * Detect a single intersection between a segment and an epipolar line.
 *  - intersection_point contains found point if found flag is set to true, its contents should be ignored otherwise
 * 	- found flag is set to false if no valid intersection is found
 */
void single_intersect_segment_line_in_detection_range(const glm::vec4 &segm,const glm::vec3 &line, const glm::vec2 &center_of_detection,const float detection_radius_squared, glm::vec2 &intersection_point, bool &found);

/**
 * Computes the epipolar line for one point, given the fundamental matrix F
 */
bool computeCorrespondEpilineSinglePoint(const glm::vec2 &starting_point_on_cur_img, const cv::Mat &F, glm::vec3 &epipolar_line, const int points_img);

vector<glm::vec3> computeCorrespondEpilineMultiplePoints(const vector<glm::vec2> &starting_points_on_cur_img, const cv::Mat &F, const int points_img);
vector<glm::vec3> computeCorrespondEpilineMultiplePoints(const vector<glm::vec2> &starting_points_on_cur_img, const cv::Mat &F, const int points_img, bool &valid_all);
/**
 * Computes minimum distance between a point and a segment.
 *
 * If the projection of the point is on the segment, distance is, as usual, the distance between
 * the point and the line the segment belongs to.
 *
 * If not, computes the distance between the point and the closest extreme of the segment.
 *
 * projection is assigned the 2D coordinates of the closest point
 */
float minimum_distance(const glm::vec4 &segm,const glm::vec2 &p, glm::vec2 &projection);

/*
* Computes minimum distance squared between a point and a segment.
*
* If the projection of the point is on the segment, distance is, as usual, the distance between
* the point and the line the segment belongs to.
*
* If not, computes the distance between the point and the closest extreme of the segment.
*
* projection is assigned the 2D coordinates of the closest point
*/
float minimum_distancesq(const glm::vec4 &segm,const glm::vec2 &p, glm::vec2 &projection);

/**
 * Computes minimum distance between a point and a segment vw.
 *
 * If the projection of the point is on the segment, distance is, as usual, the distance between
 * the point and the line the segment belongs to.
 *
 * If not, computes the distance between the point and the closest extreme of the segment.
 *
 * projection is assigned the 2D coordinates of the closest point
 */
float minimum_distance(const glm::vec2 &p,const glm::vec2 &v,const glm::vec2 &w, glm::vec2 &projection);

/*
* Computes minimum distance squared between a point and a segment vw.
*
* If the projection of the point is on the segment, distance is, as usual, the distance between
* the point and the line the segment belongs to.
*
* If not, computes the distance between the point and the closest extreme of the segment.
*
* projection is assigned the 2D coordinates of the closest point
*/
float minimum_distancesq(const glm::vec2 &p,const glm::vec2 &v,const glm::vec2 &w, glm::vec2 &projection);

bool is_visible(const glm::vec3 &pt, const CameraType &ct);

glm::vec4 compute_projection(const glm::mat4 &projection_matrix,const vec6 &segm3d);
glm::vec4 compute_projection(const glm::mat4 &projection_matrix,const glm::vec3 &segm3d_start,const glm::vec3 &segm3d_end);

glm::vec2 compute_projection(const glm::mat4 &projection_matrix,const glm::vec3 &point3d);

glm::vec2 compute_projection(const SfMData &sfmd,const int img_id,const glm::vec3 &point3d);

glm::vec2 compute_projection(const SfMData &sfmd,const int img_id,const glm::vec3 &point3d);

glm::vec3 compute_center(const glm::mat3 &rotation, const glm::vec3 &translation);

glm::vec3 compute_translation(const glm::mat3 &rotation, const glm::vec3 &center);

glm::mat4 compute_cameraMatrix(const glm::mat3 &intrinsics, const glm::mat3 &rotation, const glm::vec3 &translation);

/**
 * Compute intersection between segment and line. If an intersection is not present,
 * the function will check the endpoint on the segment directed towards the line.
 * If that point is not too far away from the line (distsq > max_close_point_distancesq)
 * that point will be selected.
 *
 * Possible results are:
 * 	- parallel
 * 		- overlapped
 * 		- not overlapped
 * 	- single intersection (only case in which intersection_found is set to true)
 * 		- close_point_selected set to true if a close point has been picked instead of an intersection
 *
 * 	Usage should first check parallel and overlapped flags, before reading the intersection field.
 */
inline void intersect_segment_line_with_close_points_on_segm(const glm::vec4 &segm,const glm::vec3 &line, const float max_close_point_distancesq, bool &parallel, bool &overlapped, bool &intersection_found, bool &close_point_selected, glm::vec2 &intersection);


/**
 * Detect a single intersection (or a close enough point) between a segment and an epipolar line.
 *  - intersection_point contains found point if found flag is set to true, its contents should be ignored otherwise
 * 	- found flag is set to false if no valid intersection is found
 */
void single_intersect_segment_line_with_close_points_on_segm_in_detection_range(const glm::vec4 &segm,const glm::vec3 &line, const glm::vec2 &center_of_detection,const float detection_radius_squared, const float max_close_point_distancesq, glm::vec2 &intersection_point, bool &found);


float distance_point_line_sq(const glm::vec2 &p2d, const glm::vec3 &line);
float distance_point_line_sq(const glm::vec3 &p3d, const vec6 &line);

float distance_point_line(const glm::vec2 &p2d, const glm::vec3 &line);
float distance_point_line(const glm::vec3 &p3d, const vec6 &line);

glm::vec3 compute_2dline(const glm::vec4 &segm);
glm::vec3 compute_2dline(const glm::vec2 &a, const glm::vec2 &b);

vec6 compute_3dline(const glm::vec3 &a, const glm::vec3 &b);

glm::vec3 convert_2Dpointonimageplane_to_3D(const glm::vec2 &image_plane_coords, const CameraType &ct);

std::vector<std::vector<std::vector<CGAL_Line3>>> compute_all_3d_lines(const SfMData &sfmd, const int starting_img_id, const int refpoint, const std::vector<SegmentEdgeManager::segment_point> &intersections, const std::vector<std::vector<std::vector<SegmentEdgeManager::segment_point>>> &correspondences);

bool cgal_3dlines_compatible(CGAL_Line3 l1, CGAL_Line3 l2);

bool cgal_3dplanes_compatible(CGAL_Plane3 p1,CGAL_Plane3 p2,CGAL_Plane3 p3);

std::vector<std::vector<CGAL_Line3>> compute_all_3d_lines_potential_match(const SfMData &sfmd, const int starting_img_id, const int refpoint, const std::vector<potential_match> &pms);

/**
 * For each plane in correspondence,
 */
vector<bool> compute_incompatible_correspondences_cgal(CGAL_Plane3 starting_plane, vector<CGAL_Plane3> correspondences_planes);

/**
 * Output format: (starting_intersections_lines, correspondence_intersections_lines)
 *
 * starting_intersections_lines : vector<CGAL_Line3>
 * 								  vector of lines resulting from intersection of starting_planes and correspondence_planes
 * 								  line is set to deprecate if intersection is not a line
 *
 * correspondence_intersections_lines : vector<vector<std::pair<int,CGAL_Line3>>>
 * 								  		vector of (correspondence_b, Line)
 * 								  		for each correspondence_a
 * 								  		-> intersection between the correspondence_a and correspondence_b
 * 								  		   if intersection is not a Line, it doesn't get added to the vector
 * 								  		NOTE: correspondence_a ignored if intersection with starting is deprecate
 * 								  		NOTE: correspondence_a < correspondence_b
 */
std::pair<vector<CGAL_Line3>, vector<vector<std::pair<int,CGAL_Line3>>>> compute_all_lines_except_incompatible(const CGAL_Plane3 &starting_plane, const std::vector<CGAL_Plane3> &corresponding_planes);

bool cgal_3dplanes_compatible(CGAL_Plane3 p1,CGAL_Plane3 p2,CGAL_Plane3 p3);
bool cgal_3dplanes_compatible(CGAL_Plane3 p1,CGAL_Plane3 p2,CGAL_Plane3 p3, CGAL_Line3 l12, CGAL_Line3 l13, CGAL_Line3 l23);

glm::vec2 middle_point(const glm::vec2 &a, const glm::vec2 &b);
glm::vec3 middle_point(const glm::vec3 &a, const glm::vec3 &b);

glm::vec2 first_plus_ratio_of_segment(const glm::vec2 &a, const glm::vec2 &b, const float ratio);

bool is_ordered_2dlinepoints(const glm::vec2 &a, const glm::vec2 &b,const glm::vec2 &c);

#endif /* INCLUDE_EDGEGRAPH3D_UTILS_GEOMETRY_GEOMETRIC_UTILITIES_HPP_ */
