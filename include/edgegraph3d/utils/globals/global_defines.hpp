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


#ifndef INCLUDE_EDGEGRAPH3D_UTILS_GLOBALS_GLOBAL_DEFINES_HPP_
#define INCLUDE_EDGEGRAPH3D_UTILS_GLOBALS_GLOBAL_DEFINES_HPP_


#define DETECTION_STARTING_RADIUS 10
#define DETECTION_CORRESPONDENCES_MULTIPLICATION_FACTOR 3
#define DETECTION_CORRESPONDENCES_MIN_RADIUS 15
#define DETECTION_CORRESPONDENCES_DETECTION_MIN_ANGLE 15
#define DETECTION_CORRESPONDENCES_DETECTION_MAX_ANGLECOS 0.965
#define DETECTION_CORRESPONDENCES_RADIUS (DETECTION_STARTING_RADIUS*DETECTION_CORRESPONDENCES_MULTIPLICATION_FACTOR)
//#define FOCUS_IMAGES_RANGE (DETECTION_CORRESPONDENCES_RADIUS*1.2)
#define FOCUS_IMAGES_RANGE 200

#define PLG3D_OUTNAME "outgraph.3dg"
#define PLGMM_OUTNAME "outplgmm.gmm"
#define PLG_STRING "optimized_plg_"
#define EDGE_COLOR Vec<unsigned char, 3>(255,255,255)

#define IMAGE_DISPLAY_TITLE "all"
//#define IMAGE_SIZE_ALL Size(720,640)
#define IMAGE_SIZE_ALL Size(1080,960)

#define COMPATIBILITY_GRAPH "potential_compatibility.graph"
#define POLYLINE_COMMUNITIES "potential_compatibility.graph_clustInfo"

#endif /* INCLUDE_EDGEGRAPH3D_UTILS_GLOBALS_GLOBAL_DEFINES_HPP_ */
