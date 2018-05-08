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


#include "geometric_utilities.hpp"

#include <boost/optional/optional.hpp>
#include <boost/variant/get.hpp>
#include <boost/variant/variant.hpp>
#include <CGAL/Kernel/global_functions_3.h>
#include <opencv2/calib3d.hpp>
#include <opencv2/core/cvdef.h>
#include <opencv2/core/mat.hpp>
#include <opencv2/core/mat.inl.hpp>
#include <opencv2/core/matx.hpp>
#include <opencv2/core/types.hpp>

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <set>

#include "SfMData.h"
#include "types_reconstructor.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "glm.hpp"


/**
 * Output: [a,b,c,d]
 * Plane is in the form: a * x + b * y + c * z + d = 0
 */
glm::vec4 plane_from_3_points_m(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3) {
	// n * (r - r0) = 0
	glm::vec3 n = glm::cross((p2-p1),(p3-p1));
	float d = -p1[0]*n[0] -p1[1]*n[1] -p1[2]*n[2];
	return glm::vec4(n[0],n[1],n[2],d);
}

std::pair<glm::vec3,float> plane_from_3_points_pairoutput_m(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3) {
	// n * (r - r0) = 0
	glm::vec3 n = glm::cross((p2-p1),(p3-p1));
	float d = -p1[0]*n[0] -p1[1]*n[1] -p1[2]*n[2];
	return std::make_pair(n,d);
}

/**
 * Output: p1, v1, v2
 * Plane is in the form: p1 + k * v1 + t * v2
 */
glm::mat3 plane_from_3_points_2vecs_form_m(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3) {
	return glm::mat3(p1,p2-p1,p3-p1);
}

/**
 * Computes the intersection between two planes:
 * 	- if parallel, set parallel true
 * 	- if overlapped (i.e. same plane), set overlapped true
 * 	- otherwise define the intersection line as:
 * 		l : x0 + k * v, k in R
 *		returned as line3d = (x0,v)
 *
 */
void planes_intersection_m(const std::pair<glm::vec3,float> &p1, const std::pair<glm::vec3,float> &p2, std::pair<glm::vec3,glm::vec3> &line3d, bool &parallel, bool &overlapped) {
	parallel = false;
	overlapped = false;

	glm::vec3 v = glm::cross(std::get<0>(p1),std::get<0>(p2)); // Direction of the intersection line

	if(v == glm::vec3(0,0,0)) {
		// Planes are parallel
		parallel = true;

		if((std::get<0>(p1)[0] / std::get<0>(p2)[0]) == (std::get<1>(p1) / std::get<1>(p2)))
			// i.e. a1 / a2 == d1 / d2
			overlapped = true;
	} else {
		// Planes are not parallel. Compute a point on the intersection to fully determine the line

		// Assuming z = 0
		glm::vec2 x0_2 = - glm::inverse(glm::mat2(std::get<0>(p1)[0],std::get<0>(p1)[1],
				std::get<0>(p2)[0],std::get<0>(p2)[1])) * glm::vec2(std::get<1>(p1),std::get<1>(p2));
		glm::vec3 x0(x0_2[0],x0_2[1],0);

		line3d = std::make_pair(x0,v);
	}
}

/**
 * Computes intersections between line described by (a,b,c) and the circle of radius r centered in (xc,yc).
 *
 * Returns vector of points of length:
 * 	0 -> no intersections
 * 	1 -> 1 intersection (i.e. line is tangent)
 * 	2 -> 2 intersections
 */
vector<glm::vec2> detect_circle_line_intersections(const glm::vec3 &line_params,const glm::vec2 &center,const float r) {
	vector<glm::vec2> result;
	float x,y,k;
	float x1,y1,x2,y2;
	float ka,kb,kc,discriminant;

	float a = line_params[0];
	float b = line_params[1];
	float c = line_params[2];

	float xc = center[0];
	float yc = center[1];

	if (b==0) {
		// Line is vertical

		x = -c/a;
		k = pow(r,2) - pow((xc + c/a),2);

		// if k < 0 no intersections are found, else
		if(k==0)
			result.push_back(glm::vec2(x,yc));
		else if (k>0) {
            y1 = yc + sqrt(k);
            y2 = yc - sqrt(k);

			result.push_back(glm::vec2(x,y1));
			result.push_back(glm::vec2(x,y2));
		}
	} else {
		// Line is not vertical


		ka=pow(a,2)+pow(b,2);
		kb=-2*pow(b,2)*xc+2*a*c+2*a*b*yc;
		kc=pow(b,2)*pow(xc,2)+pow(c,2)+2*b*c*yc+pow(b,2)*pow(yc,2)-pow(b,2)*pow(r,2);

		discriminant = pow(kb,2) - 4 * kc * ka;

		// if discriminant < 0 no intersections are found, else
		if(discriminant == 0) {
            x = - kb / (2*ka);
            y = - (a*x + c) / b;

			result.push_back(glm::vec2(x,y));
		} else if (discriminant > 0) {
            x1 = (-kb - sqrt( discriminant ))/(2*ka);
            y1 = - (a*x1 + c) / b;

            x2 = (-kb + sqrt( discriminant ))/(2*ka);
            y2 = - (a*x2 + c) / b;

            result.push_back(glm::vec2(x1,y1));
            result.push_back(glm::vec2(x2,y2));
		}
	}

	return result;
}

/**
 * Computes intersections between segment from (x1,y1) to (x2,y2) and circle of radius r centered in (xc,yc)
 */
vector<glm::vec2> detect_circle_segment_intersections(const glm::vec4 &segm,const glm::vec2 &center,const float r) {
	vector<glm::vec2> result;
	float x,y;

	float xc = center[0];
	float yc = center[1];

	float x1 = segm[0];
	float y1 = segm[1];
	float x2 = segm[2];
	float y2 = segm[3];

	// d is the vector from the segment starting point to the end
	Vec2f d(x2-x1,y2-y1);

	// f is the vector from the center of the circle to the segment start
	Vec2f f(x1-xc,y1-yc);

	float a = d.dot(d);
	float b = 2*f.dot(d) ;
	float c = f.dot(f) - r*r ;

	float discriminant = b*b-4*a*c;

	// if discriminant < 0 no intersections are found

	if (discriminant >= 0)
	{
	  // ray didn't totally miss sphere,
	  // so there is a solution to
	  // the equation.

	  discriminant = sqrt( discriminant );

	  // either solution may be on or off the ray so need to test both
	  // t1 is always the smaller value, because BOTH discriminant and
	  // a are nonnegative.
	  float t1 = (-b - discriminant)/(2*a);
	  float t2 = (-b + discriminant)/(2*a);

	  // 3x HIT cases:
	  //          -o->             --|-->  |            |  --|->
	  // Impale(t1 hit,t2 hit), Poke(t1 hit,t2>1), ExitWound(t1<0, t2 hit),

	  // 3x MISS cases:
	  //       ->  o                     o ->              | -> |
	  // FallShort (t1>1,t2>1), Past (t1<0,t2<0), CompletelyInside(t1<0, t2>1)

	  if( t1 >= 0 && t1 <= 1 )
	  {
	    // t1 is the intersection, and it's closer than t2
	    // (since t1 uses -b - discriminant)
	    // Impale, Poke
		x = x1+t1*d[0];
		y = y1+t1*d[1];
	    result.push_back(glm::vec2(x,y));
	  }

	  // here t1 didn't intersect so we are either started
	  // inside the sphere or completely past it
	  if( t2 >= 0 && t2 <= 1 )
	  {
	    // Impale,ExitWound

		x = x1+t2*d[0];
		y = y1+t2*d[1];
		result.push_back(glm::vec2(x,y));
	  }

	  // no intn: FallShort, Past, CompletelyInside
	}

	return result;
}


/**
 * Compute intersection between segment and line. Possible results are:
 * 	- parallel
 * 		- overlapped
 * 		- not overlapped
 * 	- single intersection (only case in which intersection_found is set to true)
 *
 * 	Usage should first check parallel and overlapped flags, before reading the intersection field.
 */
void intersect_segment_line(const glm::vec4 &segm,const glm::vec3 &line, bool &parallel, bool &overlapped, bool &intersection_found, glm::vec2 &intersection) {
	float num,den,t;
	glm::vec2 d(segm[2]-segm[0],segm[3]-segm[1]);

	parallel = false;
	overlapped = false;
	intersection_found = false;

	/**
	 * Segment (x1,y1)
	 * Line (a,b,c)
	 *
	 *         a*x1 + b*y1 + c
	 * t = -  _________________
	 *         a * d1 + b * d2
	 *
	 * d1=x2-x1
	 * d2=y2-y1
	 *
	 * a * d1 + b * d2 = 0 ---> same pendence
	 * a*x1 + b*y1 + + c = 0 ---> (x1,y1) belongs to the line
	 *
	 */

	num = line[0]*segm[0] + line[1]*segm[1] + line[2];
	den = line[0] * d[0] + line[1] * d[1];

	if (den!=0) {
		// a * d1 + b * d2 != 0 ---> different pendence
		t = - num/den;
		if (t >= 0 && t <= 1) {
			intersection[0] = segm[0] + t * d[0];
			intersection[1] = segm[1] + t * d[1];
			intersection_found = true;
		}
	} else {
		// a * d1 + b * d2 = 0 ---> same pendence
		parallel = true;
		overlapped = num==0;
	}
}

glm::vec3 get_line_from_ray(const ray2d &r)
{
	return glm::vec3(-r.dir[1], r.dir[0], r.dir[1] * r.start[0] - r.dir[0] * r.start[1]);
}

/**
 * Compute intersection between segment and ray. Possible results are:
 * 	- parallel
 * 		- overlapped
 * 		- not overlapped
 * 	- single intersection (only case in which intersection_found is set to true)
 *
 * 	past_center is set to true if the segment intersect the line supporting the ray before the ray starts
 *
 * 	Usage should first check parallel and overlapped flags, before reading the intersection field.
 */
void intersect_segment_ray(const glm::vec4 &segm,const ray2d &r, bool &parallel, bool &overlapped, bool &past_center, bool &intersection_found, glm::vec2 &intersection) {
	glm::vec3 line = get_line_from_ray(r);
	past_center = false;
	intersect_segment_line(segm,line, parallel, overlapped, intersection_found, intersection);
	if(intersection_found)
		if(r.dir[0] * (intersection[0] - r.start[0]) < 0) {
			intersection_found = false;
			past_center = true;
		}
}

bool aligned(const glm::vec2 &a,const glm::vec2 &b,const glm::vec2 &c)
{
	return a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y) == 0;
}

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
void intersect_segment_line_no_quasiparallel(const glm::vec4 &segm,const glm::vec3 &line, const float max_quasiparallel_angle_cos, const float max_quasiparallel_dist, bool &parallel, bool &overlapped, bool &intersection_found, bool &quasiparallel_within_distance, bool &valid, float &distance, glm::vec2 &intersection) {
	float num,den,t;
	glm::vec2 d(segm[2]-segm[0],segm[3]-segm[1]);

	valid = true;

	parallel = false;
	overlapped = false;
	intersection_found = false;
	quasiparallel_within_distance = false;

	/**
	 * Segment (x1,y1)
	 * Line (a,b,c)
	 *
	 *         a*x1 + b*y1 + c
	 * t = -  _________________
	 *         a * d1 + b * d2
	 *
	 * d1=x2-x1
	 * d2=y2-y1
	 *
	 * a * d1 + b * d2 = 0 ---> same pendence
	 * a*x1 + b*y1 + + c = 0 ---> (x1,y1) belongs to the line
	 *
	 */

	num = line[0]*segm[0] + line[1]*segm[1] + line[2];
	den = line[0] * d[0] + line[1] * d[1];

	if (den!=0) {
		// a * d1 + b * d2 != 0 ---> different pendence
		t = - num/den;
		if (t >= 0 && t <= 1) {
			intersection[0] = segm[0] + t * d[0];
			intersection[1] = segm[1] + t * d[1];
			distance = 0;
			intersection_found = true;
		}

		if(compute_anglecos(segm,line) > max_quasiparallel_angle_cos) {
			// Lines might be quasi parallel if not too distant
			if(t < 0) {
				distance = distance_point_line(glm::vec2(segm[0],segm[1]),line);
			} else if (t > 1) {
				distance = distance_point_line(glm::vec2(segm[2],segm[3]),line);
			} else {
				distance = 0;
			}
			if(distance <= max_quasiparallel_dist) {
				quasiparallel_within_distance = true;
				valid = false;
			}
		}

	} else {
		// a * d1 + b * d2 = 0 ---> same pendence
		parallel = true;
		overlapped = num==0;
		distance = distance_point_line(glm::vec2(segm[0],segm[1]),line);
		if(distance <= max_quasiparallel_dist) {
			quasiparallel_within_distance = true;
			valid = false;
		}
	}
}

bool point_in_segment_bounding_box(const glm::vec4 &segm, const glm::vec2 &pt)
{
	return ((segm[0] <= pt.x && pt.x <= segm[2]) || (segm[2] <= pt.x && pt.x <= segm[0])) && ((segm[1] <= pt.y && pt.y <= segm[3]) || (segm[3] <= pt.y && pt.y <= segm[1]));
}

void intersect_segment_segment(const glm::vec4 &segm1,const glm::vec4 &segm2, bool &parallel, bool &overlapped, bool &intersection_found, glm::vec2 &intersection)
{
	const glm::vec3 line1 = compute_2dline(segm1);
	intersect_segment_line(segm2,line1, parallel, overlapped, intersection_found, intersection);
	intersection_found = intersection_found && point_in_segment_bounding_box(segm1,intersection);
}

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
inline void intersect_segment_line_with_close_points_on_segm(const glm::vec4 &segm,const glm::vec3 &line, const float max_close_point_distancesq, bool &parallel, bool &overlapped, bool &intersection_found, bool &close_point_selected, glm::vec2 &intersection) {
	float num,den,t;
	glm::vec2 d(segm[2]-segm[0],segm[3]-segm[1]);

	parallel = false;
	overlapped = false;
	intersection_found = false;
	close_point_selected = false;

	/**
	 * Segment (x1,y1)
	 * Line (a,b,c)
	 *
	 *         a*x1 + b*y1 + c
	 * t = -  _________________
	 *         a * d1 + b * d2
	 *
	 * d1=x2-x1
	 * d2=y2-y1
	 *
	 * a * d1 + b * d2 = 0 ---> same pendence
	 * a*x1 + b*y1 + + c = 0 ---> (x1,y1) belongs to the line
	 *
	 */

	num = line[0]*segm[0] + line[1]*segm[1] + line[2];
	den = line[0] * d[0] + line[1] * d[1];

	if (den!=0) {
		// a * d1 + b * d2 != 0 ---> different pendence, not parallel. A solution may exist
		t = - num/den;
		if (t >= 0 && t <= 1) {
			intersection[0] = segm[0] + t * d[0];
			intersection[1] = segm[1] + t * d[1];
			intersection_found = true;
		} else {
			float dlensq = d[0]*d[0]+d[1]*d[1];
			// check if first segment endpoint is close enough to the line to be considered a valid close_point
			if(t<0) {
				float distdq = dlensq * t * t;

				if(dlensq < max_close_point_distancesq) {
					intersection[0] = segm[0];
					intersection[1] = segm[1];
					intersection_found = true;
					close_point_selected = true;
				}
			} else {
			// check if second segment endpoint is close enough to the line to be considered a valid close_point
				float distdq = dlensq * (t-1) * (t-1);
				if(dlensq < max_close_point_distancesq) {
					intersection[0] = segm[2];
					intersection[1] = segm[3];
					intersection_found = true;
					close_point_selected = true;
				}
			}
		}

	} else {
		// a * d1 + b * d2 = 0 ---> same pendence, parallel and potentially overlapping
		parallel = true;
		overlapped = num==0;
	}
}

/*
 * Solution:
 *
 * Parallel if a1*b2-a2*b1 = 0
 * Overlapped if parallel and b2*c1-b1*c2 = 0
 *
 * Otherwise:
 *
 *         b2*c1-b1*c2
 *  x = - -------------
 *         a1*b2-a2*b1
 *
 *         a1*c2-a2*c1
 *  y = - -------------
 *         a1*b2-a2*b1
 */
void intersect_line_line(const glm::vec3 &line_i,const glm::vec3 &line_j, bool &parallel, bool &overlapped, glm::vec2 &intersection) {
	float den = line_i[0]*line_j[1]-line_j[0]*line_i[1];
	float x_num = line_j[1]*line_i[2]-line_i[1]*line_j[2];

	parallel = den == 0;

	if(parallel)
		overlapped = x_num == 0;
	else {
		intersection[0] = - x_num / den;
		intersection[1] = - (line_i[0]*line_j[2]-line_j[0]*line_i[2]) / den;
	}
}

float squared_2d_distance(const glm::vec2 &a,const glm::vec2 &b) {
	return pow(a[0]-b[0],2) + pow(a[1]-b[1],2);
}

float squared_3d_distance(const glm::vec3 &a,const glm::vec3 &b) {
	return pow(a[0]-b[0],2) + pow(a[1]-b[1],2) + pow(a[2]-b[2],2);
}

float compute_lengthsq(const glm::vec2 &a) {
	return glm::dot(a,a);
}

float compute_lengthsq(const glm::vec3 &a) {
	return glm::dot(a,a);
}

float compute_2d_distance(const glm::vec2 &a,const glm::vec2 &b) {
	return sqrt(squared_2d_distance(a,b));
}

float compute_3d_distance(const glm::vec3 &a,const glm::vec3 &b) {
	return sqrt(squared_3d_distance(a,b));
}

float compute_anglecos_vec2_vec2(const glm::vec2 &a,const glm::vec2 &b) {
	return glm::dot(a,b) / sqrt(glm::dot(a,a)*glm::dot(b,b));
}

float compute_anglecos_vec2_vec2(const glm::vec2 &a1,const glm::vec2 &a2,const glm::vec2 &b1,const glm::vec2 &b2) {
	const glm::vec2 a = a2-a1;
	const glm::vec2 b = b2-b1;

	return compute_anglecos_vec2_vec2(a,b);
}

/**
 * Detect a single intersection between a segment and an epipolar line.
 *  - intersection_point contains found point if found flag is set to true, its contents should be ignored otherwise
 * 	- found flag is set to false if no valid intersection is found
 */
void single_intersect_segment_line_in_detection_range(const glm::vec4 &segm,const glm::vec3 &line, const glm::vec2 &center_of_detection,const float detection_radius_squared, glm::vec2 &intersection_point, bool &found) {
	bool parallel,overlapped,intersection_found;

	intersect_segment_line(segm,line, parallel, overlapped, intersection_found, intersection_point);

	found = intersection_found && (squared_2d_distance(center_of_detection,intersection_point) < detection_radius_squared);
}

glm::vec2 get_2d_direction(const glm::vec4 &segm)
{
	return glm::vec2(segm[2]-segm[0],segm[3]-segm[1]);
}

glm::vec2 get_2d_direction(const glm::vec3 &line)
{
	if(line[1] == 0)
		return glm::vec2(0.0,1.0);
	else
		return glm::vec2(1.0, -line[0]/line[1]);
}

float compute_anglecos(const glm::vec4 &segm,const glm::vec3 &line) {
	return compute_anglecos_vec2_vec2(get_2d_direction(segm),get_2d_direction(line));
}

/**
 * Detect a single intersection (or a close enough point) between a segment and an epipolar line.
 *  - intersection_point contains found point if found flag is set to true, its contents should be ignored otherwise
 * 	- found flag is set to false if no valid intersection is found
 */
void single_intersect_segment_line_with_close_points_on_segm_in_detection_range(const glm::vec4 &segm,const glm::vec3 &line, const glm::vec2 &center_of_detection,const float detection_radius_squared, const float max_close_point_distancesq, glm::vec2 &intersection_point, bool &found) {
	bool parallel,overlapped,intersection_found;
	bool close_point_selected;

	intersect_segment_line_with_close_points_on_segm(segm,line, max_close_point_distancesq, parallel, overlapped, intersection_found, close_point_selected, intersection_point);

	found = intersection_found && (squared_2d_distance(center_of_detection,intersection_point) < detection_radius_squared);
}

void t_relativeCameraMotion(const Mat &R1,
                          const Mat &t1,
                          const Mat &R2,
                          const Mat &t2,
						  Mat &R,
						  Mat &t) {
  R = R2 * R1.t();
  t = t2 - R * t1;
}

template<typename T>
void t_essentialFromRt(const Mat &R1,
        const Mat &t1,
        const Mat &R2,
        const Mat &t2,
		  Mat &E) {
	Mat R(3, 3, CV_32FC2);
	Mat t(3, 1, CV_32FC2);
	//t_relativeCameraMotion(R1, t1, R2, t2, R, t);
	  R = R2 * R1.t();
	  t = t2 - R * t1;
  Mat cross_prod_t(3, 3, CV_32FC2);

  cross_prod_t.at<float>(0,0)=0;
  cross_prod_t.at<float>(0,1)=-(t.at<float>(2,0));
  cross_prod_t.at<float>(0,2)=t.at<float>(1,0);

  cross_prod_t.at<float>(1,0)=t.at<float>(2,0);
  cross_prod_t.at<float>(1,1)=0;
  cross_prod_t.at<float>(1,2)=-(t.at<float>(0,0));

  cross_prod_t.at<float>(2,0)=-(t.at<float>(1,0));
  cross_prod_t.at<float>(2,1)=t.at<float>(0,0);
  cross_prod_t.at<float>(2,2)=0;

	cout << "R\n";
	print_Mat<float>(R);

	cout << "T\n";
	print_Mat<float>(t);

	cout << "Tx\n";
	print_Mat<float>(cross_prod_t);

  E = cross_prod_t * R;

  cross_prod_t.release();
}

void findFundamentalMatrixFromRt(const glm::vec3 &t1,const glm::mat3 &r1,glm::mat3 i1, glm::vec3 t2,glm::mat3 r2, glm::mat3 i2, cv::Mat &F) {
	glm::mat3 E_,F_;

	glm::vec3 t;
	glm::mat3 R,Tx;


	  R = r2 * glm::transpose(r1);
	  t = t2 - R * t1;

	  Tx[0][0]=0;
	  Tx[0][1]=-t[2];
	  Tx[0][2]=t[1];

	  Tx[1][0]=t[2];
	  Tx[1][1]=0;
	  Tx[1][2]=-t[0];

	  Tx[2][0]=-t[1];
	  Tx[2][1]=t[0];
	  Tx[2][2]=0;

	  E_ = Tx * R;

	  F_ = glm::transpose(glm::inverse(i1)) * E_ * glm::inverse(i2);

	  convert_from_glm_mat3_to_cv_Mat3f(F_,F);
}

void findFundamentalMatrixFromRt2(const glm::vec3 &t1,const glm::mat3 &r1,glm::mat3 i1, glm::vec3 t2,glm::mat3 r2, glm::mat3 i2, cv::Mat &F) {

	cv::Mat mt1(3, 1, CV_32FC2);
	cv::Mat mr1(3, 3, CV_32FC2);
	cv::Mat mi1(3, 3, CV_32FC2);
	cv::Mat mt2(3, 1, CV_32FC2);
	cv::Mat mr2(3, 3, CV_32FC2);
	cv::Mat mi2(3, 3, CV_32FC2);

	convert_from_glm_vec3_to_cv_Mat(t1,mt1);
	convert_from_glm_mat3_to_cv_Mat3f(r1,mr1);
	convert_from_glm_mat3_to_cv_Mat3f(i1,mi1);
	convert_from_glm_vec3_to_cv_Mat(t2,mt2);
	convert_from_glm_mat3_to_cv_Mat3f(r2,mr2);
	convert_from_glm_mat3_to_cv_Mat3f(i2,mi2);

	cv::Mat_<float> E(3, 3, CV_32FC2);

	E.at<float>(1,1)=3.2;

	cout << "E\n";
	print_Mat<float>(E);

	t_essentialFromRt<float>(mr1,
				mt1,
				mr2,
				mt2,E);

	cout << "E\n";
	print_Mat<float>(E);


	F = mi1.inv().t() * E * mi2.inv();

	cout << "F\n";
	print_Mat<float>(F);
}

#define FAKE_CORRESPONDENCES_AMOUNT 1000

#define MIN_CORRESPONDENCES_AMOUNT 10

void findFundamentalMatrixFromPoints(const vector<set<int>> &points_on_images,const SfMData &sfmd, int i, int j, cv::Mat &F) {
	set<int> points_on_both = find_points_on_both_images(points_on_images, i, j);

	if(points_on_both.size()>=MIN_CORRESPONDENCES_AMOUNT) {
		vector<Point2f> points_i;
		vector<Point2f> points_j;
		if(points_on_both.size()>0) {
			glm::vec2 p2d;
			bool found;
			for(const auto common_point_id : points_on_both) {
				get_2d_coordinates_of_point_on_image(sfmd, i, common_point_id, p2d, found);
				points_i.push_back(Point2f(p2d[0],p2d[1]));

				get_2d_coordinates_of_point_on_image(sfmd, j, common_point_id, p2d, found);
				points_j.push_back(Point2f(p2d[0],p2d[1]));
			}
		} else {
			vector<pair<glm::vec2,glm::vec2>> fake_correspondences = find_fakepoints_on_both_images(sfmd, i, j, FAKE_CORRESPONDENCES_AMOUNT);
			for(const auto &pp:fake_correspondences) {
				points_i.push_back(Point2f(pp.first[0],pp.first[1]));
				points_j.push_back(Point2f(pp.second[0],pp.second[1]));
			}
		}

		F=findFundamentalMat(points_i, points_j,FM_LMEDS);
	} else
		F=Mat(1,1,CV_32F);
}


void generate_fundamental_matrix_from_Rt(const SfMData &sfmd, int i, int j, cv::Mat &F) {
	findFundamentalMatrixFromRt(sfmd.camerasList_[i].translation,
			sfmd.camerasList_[i].rotation,
			sfmd.camerasList_[i].intrinsics,
			sfmd.camerasList_[j].translation,
			sfmd.camerasList_[j].rotation,
			sfmd.camerasList_[j].intrinsics, F);
}

cv::Mat** generate_all_fundamental_matrices_from_Rt(const SfMData &sfmd) {
	cv::Mat** all_fms = create_2D_Mat_array((unsigned long) sfmd.numCameras_,(unsigned long) sfmd.numCameras_);

	for(int i=0; i < sfmd.numCameras_;i++)
		for(int j=0;j < sfmd.numCameras_;j++)
			if(i!=j)
				generate_fundamental_matrix_from_Rt(sfmd, i, j, all_fms[i][j]);

	return all_fms;
}


cv::Mat** generate_all_fundamental_matrices_from_Points(const SfMData &sfmd) {
	vector<set<int>> points_on_images = get_point_sets_on_images(sfmd);

	cv::Mat** all_fms = create_2D_Mat_array((unsigned long) sfmd.numCameras_,(unsigned long) sfmd.numCameras_);

	for(int i=0; i < sfmd.numCameras_;i++)
		for(int j=0;j < sfmd.numCameras_;j++)
			if(i!=j)
				findFundamentalMatrixFromPoints(points_on_images,sfmd, i, j, all_fms[i][j]);

	return all_fms;
}

cv::Mat** generate_all_fundamental_matrices(const SfMData &sfmd) {
	return generate_all_fundamental_matrices_from_Points(sfmd);
}



bool computeCorrespondEpilineSinglePoint(const glm::vec2 &starting_point_on_cur_img, const cv::Mat &F, glm::vec3 &epipolar_line, const int points_img) {
	bool valid =true;
	if(F.rows == 3 && F.cols == 3) {
		try {
			vector<Point2f> points;
			points.push_back(Point2f(starting_point_on_cur_img[0],starting_point_on_cur_img[1]));

			vector<Vec3f> eplines;
			cv::computeCorrespondEpilines(cv::Mat(points), points_img, F, eplines);

			epipolar_line[0] = eplines[0][0];
			epipolar_line[1] = eplines[0][1];
			epipolar_line[2] = eplines[0][2];
		} catch(exception& e) {
			valid=false;
		}
	} else
		valid=false;
	return valid;
}

vector<glm::vec3> computeCorrespondEpilineMultiplePoints(const vector<glm::vec2> &starting_points_on_cur_img, const cv::Mat &F, const int points_img) {
	vector<glm::vec3> epilines;
	glm::vec3 epiline;
	for(const auto &starting_point_on_cur_img : starting_points_on_cur_img) {
		if(computeCorrespondEpilineSinglePoint(starting_point_on_cur_img, F, epiline, points_img))
			epilines.push_back(epiline);
	}
	return epilines;
}

vector<glm::vec3> computeCorrespondEpilineMultiplePoints(const vector<glm::vec2> &starting_points_on_cur_img, const cv::Mat &F, const int points_img, bool &valid_all) {
	vector<glm::vec3> epilines;
	valid_all=true;
	glm::vec3 epiline;
	for(const auto &starting_point_on_cur_img : starting_points_on_cur_img) {
		if(computeCorrespondEpilineSinglePoint(starting_point_on_cur_img, F, epiline, points_img))
			epilines.push_back(epiline);
		else
			valid_all=false;
	}
	return epilines;
}


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
float minimum_distance(const glm::vec4 &segm,const glm::vec2 &p, glm::vec2 &projection) {
	return sqrt(minimum_distancesq(segm, p, projection));
}

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
float minimum_distancesq(const glm::vec4 &segm,const glm::vec2 &p, glm::vec2 &projection) {
	glm::vec2 v,w;

	v[0]=segm[0];
	v[1]=segm[1];
	w[0]=segm[2];
	w[1]=segm[3];

  // Return minimum distance between line segment vw and point p
  const float l2 = squared_2d_distance(v, w);  // i.e. |w-v|^2 -  avoid a sqrt
  if (l2 == 0.0) {
	  projection = v;
	  return squared_2d_distance(p, v);   // v == w case
  }
  // Consider the line extending the segment, parameterized as v + t (w - v).
  // We find projection of point p onto the line.
  // It falls where t = [(p-v) . (w-v)] / |w-v|^2
  // We clamp t from [0,1] to handle points outside the segment vw.
  const float t = max<float>(0, min<float>(1, dot(p - v, w - v) / l2));
  projection = v + t * (w - v);  // Projection falls on the segment
  return squared_2d_distance(p, projection);
}

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
float minimum_distance(const glm::vec2 &p,const glm::vec2 &v,const glm::vec2 &w, glm::vec2 &projection) {
	return sqrt(minimum_distancesq(p, v, w, projection));
}

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
float minimum_distancesq(const glm::vec2 &p,const glm::vec2 &v,const glm::vec2 &w, glm::vec2 &projection) {
  // Return minimum distance between line segment vw and point p
  const float l2 = squared_2d_distance(v, w);  // i.e. |w-v|^2 -  avoid a sqrt
  if (l2 == 0.0) {
	  projection = v;
	  return squared_2d_distance(p, v);   // v == w case
  }
  // Consider the line extending the segment, parameterized as v + t (w - v).
  // We find projection of point p onto the line.
  // It falls where t = [(p-v) . (w-v)] / |w-v|^2
  // We clamp t from [0,1] to handle points outside the segment vw.
  const float t = max<float>(0, min<float>(1, dot(p - v, w - v) / l2));
  projection = v + t * (w - v);  // Projection falls on the segment
  return squared_2d_distance(p, projection);
}

bool is_visible(const glm::vec3 &pt, const CameraType &ct) {
	glm::vec2 pt2d = compute_projection(ct.cameraMatrix,pt);
	return pt2d.x >= 0 && pt2d.x < ct.imageWidth && pt2d.y >= 0 && pt2d.y < ct.imageHeight;
}

glm::vec4 compute_projection(const glm::mat4 &projection_matrix,const vec6 &segm3d) {
	glm::vec2 start = compute_projection(projection_matrix,segm3d.start());
	glm::vec2 end = compute_projection(projection_matrix,segm3d.end());
	return glm::vec4(start,end);
}

glm::vec4 compute_projection(const glm::mat4 &projection_matrix,const glm::vec3 &segm3d_start,const glm::vec3 &segm3d_end) {
	glm::vec2 start = compute_projection(projection_matrix,segm3d_start);
	glm::vec2 end = compute_projection(projection_matrix,segm3d_end);
	return glm::vec4(start,end);
}

glm::vec2 compute_projection(const glm::mat4 &projection_matrix,const glm::vec3 &point3d) {
	glm::vec4 point4d = glm::vec4(point3d[0],point3d[1],point3d[2],1.0) * projection_matrix;

	return glm::vec2(point4d[0]/point4d[2],point4d[1]/point4d[2]);
}

glm::vec2 compute_projection(const SfMData &sfmd,const int img_id,const glm::vec3 &point3d) {
	const glm::mat4 &projection_matrix = sfmd.camerasList_[img_id].cameraMatrix;

	return compute_projection(sfmd.camerasList_[img_id].cameraMatrix, point3d);
}

glm::vec3 compute_center(const glm::mat3 &rotation, const glm::vec3 &translation) {
	return - translation * glm::transpose(rotation);
}

glm::vec3 compute_translation(const glm::mat3 &rotation, const glm::vec3 &center) {
	return -rotation * center;
}

glm::mat4 compute_cameraMatrix(const glm::mat3 &intrinsics, const glm::mat3 &rotation, const glm::vec3 &translation) {
	return get_rt4x4(rotation,translation) * glm::mat4(intrinsics[0][0],intrinsics[0][1],intrinsics[0][2],0,intrinsics[1][0],intrinsics[1][1],intrinsics[1][2],0,intrinsics[2][0],intrinsics[2][1],intrinsics[2][2],0,0,0,0,0);
}

float distance_point_line_sq(const glm::vec2 &p2d, const glm::vec3 &line) {
	float den = line[0] * p2d[0] + line[1] * p2d[1] + line[2];
	den*=den;
	return den/(line[0] * line[0] + line[1] * line[1]);
}

float distance_point_line_sq(const glm::vec3 &p3d, const vec6 &line) {
	return compute_lengthsq(glm::cross(line.direction(),(line.start()-p3d))) / compute_lengthsq(line.direction());
}

float distance_point_line(const glm::vec2 &p2d, const glm::vec3 &line) {
	return sqrt(distance_point_line_sq(p2d, line));
}

float distance_point_line(const glm::vec3 &p3d, const vec6 &line) {
	return sqrt(distance_point_line_sq(p3d, line));
}


/**
 * Compute matrix composed by first column of R, second column of R, and T
 */
glm::mat3 get_r1r2t(const glm::mat3 &r, const glm::vec3 &t) {
	return glm::mat3(
			r[0][0],r[0][1],t[0],
			r[1][0],r[1][1],t[1],
			r[2][0],r[2][1],t[2]);
}

glm::mat3x4 get_rt(const glm::mat3 &r, const glm::vec3 &t) {
	return glm::mat3x4(
			r[0][0],r[0][1],r[0][2],t[0],
			r[1][0],r[1][1],r[1][2],t[1],
			r[2][0],r[2][1],r[2][2],t[2]);
}

glm::mat3 get_H(const CameraType &ct) {
	return ct.intrinsics * get_r1r2t(ct.rotation,ct.translation) / ct.translation[2];
	//return ct.intrinsics * get_r1r2t(ct.rotation,ct.translation);
}

glm::vec3 convert_2Dpointonimageplane_to_3D(const glm::vec2 &image_plane_coords, const CameraType &ct) {

	const glm::mat3 &rotation = ct.rotation;
	const glm::vec3 &translation = ct.translation;

	cout << "Intrinsics: \n" << ct.intrinsics << "\n";
	cout << "R: \n" << ct.rotation << "\n";
	cout << "t: \n" << ct.translation << "\n";
	cout << "CM: \n" << ct.cameraMatrix << "\n";

	cout << "\n\n2D point: " << image_plane_coords << "\n";

	double u = image_plane_coords[0];
	double v = image_plane_coords[1];

	float f = ct .intrinsics[0][0];
	float uc = ct .intrinsics[0][2];
	float vc = ct .intrinsics[1][2];

	glm::vec4 cameracoords(u-uc,v-vc,f,1);

	cout << "cameracoords: \n" << cameracoords << "\n";

	cout << "Rt: \n" << get_rt4x4(rotation,translation) << "\n";

	glm::mat4 invrt = glm::inverse(get_rt4x4(rotation,translation));
	cout << "Rt^-1: \n" << invrt << "\n";

	glm::vec4 pt3DHomo = cameracoords * invrt;
	glm::vec3 p3d = glm::vec3(pt3DHomo[0] / pt3DHomo[3], pt3DHomo[1] / pt3DHomo[3], pt3DHomo[2] / pt3DHomo[3]);

	glm::vec2 prj_of_3dpoint = compute_projection(ct.cameraMatrix,p3d);
	cout << "Projection of computed 3D point: (" << prj_of_3dpoint[0] << "," << prj_of_3dpoint[1] << "), originated from image plane point: (" << image_plane_coords[0] << "," << image_plane_coords[1] << ")" << endl;

	return p3d;
}

std::pair<glm::vec3,glm::vec3> convert_2Dsegmentonimageplane_to_3D(const SfMData &sfmd, const int img_id, const glm::vec4 &segm2d) {
	glm::vec2 extreme2D_A(segm2d[0],segm2d[1]);
	glm::vec3 extreme3D_A = convert_2Dpointonimageplane_to_3D(extreme2D_A,sfmd.camerasList_[img_id]);

	glm::vec2 extreme2D_B(segm2d[2],segm2d[3]);
	glm::vec3 extreme3D_B = convert_2Dpointonimageplane_to_3D(extreme2D_B,sfmd.camerasList_[img_id]);

	return make_pair(extreme3D_A,extreme3D_B);
}

CGAL_Plane3 cgal_plane_from_glmvec3points(const glm::vec3 &point_A,const glm::vec3 &point_B,const glm::vec3 &point_C) {
	return CGAL_Plane3(CGAL_Point3(point_A[0],point_A[1],point_A[2]),CGAL_Point3(point_B[0],point_B[1],point_B[2]),CGAL_Point3(point_C[0],point_C[1],point_C[2]));
}

CGAL_Plane3 compute_plane_from_2d_edge(const SfMData &sfmd, const int img_id, const glm::vec4 &segm2d) {
	glm::vec3 cam_coords = sfmd.camerasList_[img_id].center;
	std::pair<glm::vec3,glm::vec3> segment_coords = convert_2Dsegmentonimageplane_to_3D(sfmd, img_id, segm2d);
	return cgal_plane_from_glmvec3points(cam_coords,segment_coords.first, segment_coords.second);
}

/**
 * Returns regular plane intersection if its a line. If the planes are coincident, or don't intersect, return degenerate line
 */
CGAL_Line3 planes_intersection_line_or_degenerate(const CGAL_Plane3 &p1, const CGAL_Plane3 &p2) {
	boost::optional<boost::variant<CGAL_Line3,CGAL_Plane3>> intersection_result = CGAL::intersection(p1,p2);
	if(intersection_result) {
		const CGAL_Line3* l = boost::get<CGAL_Line3>(&*intersection_result);
		if(l)
			return *l;
	}

	// Put a degenerate line in result
	return CGAL_Line3(CGAL_Point3(0,0,0),CGAL_Point3(0,0,0));
}

std::vector<std::vector<CGAL_Line3>> compute_3d_lines(const SfMData &sfmd, const int starting_img_id, const int refpoint, const SegmentEdgeManager::segment_point &starting_intersection, const std::vector<std::vector<SegmentEdgeManager::segment_point>> &correspondences) {
	CGAL_Plane3 starting_edge_plane = compute_plane_from_2d_edge(sfmd,starting_img_id,starting_intersection.segment_extremes);
	std::vector<std::vector<CGAL_Line3>> res;

	for(int i=0; i < sfmd.camViewingPointN_[refpoint].size(); i++) {
		const int img_id = sfmd.camViewingPointN_[refpoint][i];
		std::vector<CGAL_Line3> current_img_res;

		if(img_id != starting_img_id) {
			const std::vector<SegmentEdgeManager::segment_point> &current_correspondences = correspondences[img_id];
			for(const auto &correspondence : current_correspondences) {
				CGAL_Plane3 current_edge_plane = compute_plane_from_2d_edge(sfmd,img_id,correspondence.segment_extremes);

				current_img_res.push_back(planes_intersection_line_or_degenerate(starting_edge_plane,current_edge_plane));
			}
		} else {
			// Put a degenerate line in result
			current_img_res.push_back(CGAL_Line3(CGAL_Point3(0,0,0),CGAL_Point3(0,0,0)));
		}

		res.push_back(current_img_res);
	}

	return res;
}

std::vector<CGAL_Line3> compute_3d_lines_potential_match(const SfMData &sfmd, const int starting_img_id, const int refpoint, const potential_match &pm) {
	const SegmentEdgeManager::segment_point &starting_intersection = pm.starting_intersection;

	CGAL_Plane3 starting_edge_plane = compute_plane_from_2d_edge(sfmd,starting_img_id,starting_intersection.segment_extremes);
	std::vector<CGAL_Line3> res;

	for(int i=0; i < sfmd.camViewingPointN_[refpoint].size(); i++) {
		const int img_id = sfmd.camViewingPointN_[refpoint][i];

		if(img_id != starting_img_id) {
			const std::vector<uint> &current_correspondences_ids = pm.correspondences[i];
			for(const auto current_correspondence_id : current_correspondences_ids) {
				const SegmentEdgeManager::segment_point &current_correspondence = pm.correspondences_data[current_correspondence_id];
				CGAL_Plane3 current_edge_plane = compute_plane_from_2d_edge(sfmd,img_id,current_correspondence.segment_extremes);

				// Compute planes intersection
				res.push_back(planes_intersection_line_or_degenerate(starting_edge_plane,current_edge_plane));
			}
		} else {
			const std::vector<uint> &current_correspondences_ids = pm.correspondences[i];
			if(current_correspondences_ids.size() > 0)
				// Put a degenerate line in result
				res.push_back(CGAL_Line3(CGAL_Point3(0,0,0),CGAL_Point3(0,0,0)));
		}
	}

	return res;
}

std::vector<std::vector<std::vector<CGAL_Line3>>> compute_all_3d_lines(const SfMData &sfmd, const int starting_img_id, const int refpoint, const std::vector<SegmentEdgeManager::segment_point> &intersections, const std::vector<std::vector<std::vector<SegmentEdgeManager::segment_point>>> &correspondences) {
	std::vector<std::vector<std::vector<CGAL_Line3>>> res;

	for(int intersection_index=0; intersection_index < intersections.size();intersection_index++) {
		const SegmentEdgeManager::segment_point &starting_intersection = intersections[intersection_index];
		const std::vector<std::vector<SegmentEdgeManager::segment_point>> &current_correspondences = correspondences[intersection_index];

		res.push_back(compute_3d_lines(sfmd,starting_img_id, refpoint,starting_intersection,current_correspondences));
	}

	return res;
}

std::vector<std::vector<CGAL_Line3>> compute_all_3d_lines_potential_match(const SfMData &sfmd, const int starting_img_id, const int refpoint, const std::vector<potential_match> &pms) {
	std::vector<std::vector<CGAL_Line3>> res;

	for(const auto &pm: pms) {
		res.push_back(compute_3d_lines_potential_match(sfmd,starting_img_id, refpoint,pm));
	}

	return res;
}

inline bool cgal_3dlines_compatible_3d_angle(CGAL_Line3 l1, CGAL_Line3 l2) {
	if(l1.is_degenerate() || l2.is_degenerate())
		return false;

	CGAL_Vector3 v1 = l1.to_vector();
	CGAL_Vector3 v2 = l2.to_vector();

	double angle_cos = v1 * v2 / sqrt(v1.squared_length() * v2.squared_length());

	cout << "Computed 3D angle cos: " << angle_cos << endl;

	return abs(angle_cos) <= MAX_LINE3DANGLE_COS;
}

/**
 * output (starting_segment_plane, correspondences_planes)
 */
std::pair<CGAL_Plane3, std::vector<CGAL_Plane3>> compute_3d_planes_potential_match(const SfMData &sfmd, const int starting_img_id, const int refpoint, const potential_match &pm) {
	const SegmentEdgeManager::segment_point &starting_intersection = pm.starting_intersection;

	CGAL_Plane3 starting_edge_plane = compute_plane_from_2d_edge(sfmd,starting_img_id,starting_intersection.segment_extremes);

	std::vector<CGAL_Plane3> res;

	for(int i=0; i < sfmd.camViewingPointN_[refpoint].size(); i++) {
		const int img_id = sfmd.camViewingPointN_[refpoint][i];

		if(img_id != starting_img_id) {
			const std::vector<uint> &current_correspondences_ids = pm.correspondences[i];
			for(const auto current_correspondence_id : current_correspondences_ids) {
				const SegmentEdgeManager::segment_point &current_correspondence = pm.correspondences_data[current_correspondence_id];
				CGAL_Plane3 current_edge_plane = compute_plane_from_2d_edge(sfmd,img_id,current_correspondence.segment_extremes);

				res.push_back(current_edge_plane);
			}
		} else {
			const std::vector<uint> &current_correspondences_ids = pm.correspondences[i];
			if(current_correspondences_ids.size() > 0)
				// Put a degenerate line in result
				res.push_back(starting_edge_plane);
		}
	}

	return make_pair(starting_edge_plane,res);
}

#define ALMOST_PERPENDICULAR_MAXCOS 0.1

inline bool cgal_vectors_almost_perpendicular(CGAL_Vector3 v1, CGAL_Vector3 v2) {
	return abs(v1 * v2 / sqrt(v1.squared_length() * v2.squared_length())) < ALMOST_PERPENDICULAR_MAXCOS;
}

bool cgal_3dplanes_compatible(CGAL_Plane3 p1,CGAL_Plane3 p2,CGAL_Plane3 p3) {
	CGAL_Line3 l12 = planes_intersection_line_or_degenerate(p1,p2);
	if(l12.is_degenerate())
		return false;
	CGAL_Line3 l13 = planes_intersection_line_or_degenerate(p1,p3);
	if(l13.is_degenerate())
		return false;
	CGAL_Line3 l23 = planes_intersection_line_or_degenerate(p2,p3);
	if(l23.is_degenerate())
		return false;

	CGAL_Vector3 v12 = l12.to_vector();
	//v12 /= sqrt(v12.squared_length());
	CGAL_Vector3 v13 = l13.to_vector();
	//v13 /= sqrt(v13.squared_length());
	CGAL_Vector3 v23 = l23.to_vector();
	//v23 /= sqrt(v23.squared_length());

	CGAL_Vector3 v1 = p1.orthogonal_vector();
	//v1 /= v1.squared_length();
	CGAL_Vector3 v2 = p2.orthogonal_vector();
	//v2 /= v2.squared_length();
	CGAL_Vector3 v3 = p3.orthogonal_vector();
	//v3 /= v3.squared_length();

	// Plane normals and lines should be ~ perpendicular for the lines to belong to the plane

	return cgal_vectors_almost_perpendicular(v1,v23) && cgal_vectors_almost_perpendicular(v2,v13) && cgal_vectors_almost_perpendicular(v3,v12);

}

bool cgal_3dplanes_compatible(CGAL_Plane3 p1,CGAL_Plane3 p2,CGAL_Plane3 p3, CGAL_Line3 l12, CGAL_Line3 l13, CGAL_Line3 l23) {
	CGAL_Vector3 v12 = l12.to_vector();
	//v12 /= sqrt(v12.squared_length());
	CGAL_Vector3 v13 = l13.to_vector();
	//v13 /= sqrt(v13.squared_length());
	CGAL_Vector3 v23 = l23.to_vector();
	//v23 /= sqrt(v23.squared_length());

	CGAL_Vector3 v1 = p1.orthogonal_vector();
	//v1 /= v1.squared_length();
	CGAL_Vector3 v2 = p2.orthogonal_vector();
	//v2 /= v2.squared_length();
	CGAL_Vector3 v3 = p3.orthogonal_vector();
	//v3 /= v3.squared_length();

	// Plane normals and lines should be ~ perpendicular for the lines to belong to the plane

	return cgal_vectors_almost_perpendicular(v1,v23) && cgal_vectors_almost_perpendicular(v2,v13) && cgal_vectors_almost_perpendicular(v3,v12);

}

bool cgal_3dlines_compatible(CGAL_Line3 l1, CGAL_Line3 l2) {
	return cgal_3dlines_compatible_3d_angle(l1,l2);
}

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
std::pair<vector<CGAL_Line3>, vector<vector<std::pair<int,CGAL_Line3>>>> compute_all_lines_except_incompatible(const CGAL_Plane3 &starting_plane, const std::vector<CGAL_Plane3> &corresponding_planes) {
	const int amount_of_correspondences = corresponding_planes.size();

	vector<CGAL_Line3> s_intersection_lines;
	vector<vector<std::pair<int,CGAL_Line3>>> c_intersection_lines;

	for(int i=0; i < amount_of_correspondences; i++)
		s_intersection_lines.push_back(planes_intersection_line_or_degenerate(starting_plane,corresponding_planes[i]));

	for(int correspondence_a=0; correspondence_a < amount_of_correspondences; correspondence_a++) {
		vector<std::pair<int,CGAL_Line3>> correspondence_a_compatibles;

		if(!(s_intersection_lines[correspondence_a].is_degenerate()))
			for(int correspondence_b=correspondence_a+1; correspondence_b < amount_of_correspondences; correspondence_b++)
				if(!(s_intersection_lines[correspondence_b].is_degenerate())) {
					CGAL_Line3 intersection_ab = planes_intersection_line_or_degenerate(corresponding_planes[correspondence_a],corresponding_planes[correspondence_b]);
					if(!(intersection_ab.is_degenerate()))
						correspondence_a_compatibles.push_back(make_pair(correspondence_b, intersection_ab));
				}

		c_intersection_lines.push_back(correspondence_a_compatibles);
	}

	return make_pair(s_intersection_lines,c_intersection_lines);
}

glm::vec3 compute_2dline(const glm::vec4 &segm) {
	return compute_2dline(glm::vec2(segm[0],segm[1]),glm::vec2(segm[2],segm[3]));
}

glm::vec3 compute_2dline(const glm::vec2 &a, const glm::vec2 &b) {
	glm::vec3 res;

	if(a.x == b.x) {
		// vertical line
		res = glm::vec3(1,0,-a.x);
	} else {
		float m = (b.y - a.y) / (b.x - a.x);
		float q = a.y - m * a.x;
		res = glm::vec3(m,-1.0,q);
	}

	return res;
}

vec6 compute_3dline(const glm::vec3 &a, const glm::vec3 &b) {
	return vec6(a[0],a[1],a[2],b[0],b[1],b[2]);
}



glm::vec2 middle_point(const glm::vec2 &a, const glm::vec2 &b) {
	return glm::vec2((a[0]+b[0])/2,(a[1]+b[1])/2);
}

glm::vec3 middle_point(const glm::vec3 &a, const glm::vec3 &b) {
	return glm::vec3((a[0]+b[0])/2,(a[1]+b[1])/2,(a[2]+b[2])/2);
}

glm::vec2 first_plus_ratio_of_segment(const glm::vec2 &a, const glm::vec2 &b, const float ratio) {
	return glm::vec2(a[0]+ratio * (b[0]-a[0]), a[1]+ratio * (b[1]-a[1]));
}

// assuming a,b,c are on a 2D line, returns true if a <= b <= c or c <= b <= a
bool is_ordered_2dlinepoints(const glm::vec2 &a, const glm::vec2 &b,const glm::vec2 &c) {
	return (b.x-a.x)*(c.x-b.x) > 0 || (b.y-a.y)*(c.y-b.y) > 0 || (a==b || b==c);
}
