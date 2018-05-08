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


#ifndef INCLUDE_EDGEGRAPH3D_IO_OUTPUT_OUTPUT_UTILITIES_HPP_
#define INCLUDE_EDGEGRAPH3D_IO_OUTPUT_OUTPUT_UTILITIES_HPP_

#include <opencv2/core/mat.hpp>
#include <iostream>
#include <string>
#include <vector>

#include "polyline_graph_2d.hpp"
#include "polyline_graph_2d_hmap_impl.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"
#include "glm.hpp"


void print_new_point_to_file(const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &p3d, ofstream &outfile, const ulong current_point_id, const ulong last_point_id);

void print_new_points_to_file(const SfMData &sfmd, const vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &p3ds, const string &outfilename);

void add_3dpoints_to_sfmd(SfMData &sfmd, const vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &p3ds);

void write_sfmd_plus_points(const char *original_sfm_data_file, vector<Mat> &imgs, const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &p3ds, const char *out_folder, const char *newSfMname);

#endif /* INCLUDE_EDGEGRAPH3D_IO_OUTPUT_OUTPUT_UTILITIES_HPP_ */
