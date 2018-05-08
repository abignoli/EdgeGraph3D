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


#ifndef INCLUDE_EDGEGRAPH3D_PLGS_PLG_HANDLING_HPP_
#define INCLUDE_EDGEGRAPH3D_PLGS_PLG_HANDLING_HPP_

#include <opencv2/core/cvstd.hpp>
#include <opencv2/core/mat.hpp>
#include <opencv2/core/matx.hpp>
#include <vector>

#include "polyline_graph_2d_hmap_impl.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"

void write_plgs(const vector<PolyLineGraph2DHMapImpl> &plgs, const SfMData &sfmd, char *em_out_folder,  const String &s);

vector<PolyLineGraph2DHMapImpl> compute_and_write_plgs(const vector<Mat> &edge_imgs, const Vec<unsigned char, 3> &edge_color, const SfMData &sfmd, char *em_out_folder,  const String &s);

vector<PolyLineGraph2DHMapImpl> read_plgs(const SfMData &sfmd, char *em_out_folder,  const String &s);

void transform_plgs(vector<PolyLineGraph2DHMapImpl> &plgs);

#endif /* INCLUDE_EDGEGRAPH3D_PLGS_PLG_HANDLING_HPP_ */
