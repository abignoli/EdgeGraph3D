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


#ifndef INCLUDE_EDGEGRAPH3D_IO_OUTPUT_OUTPUT_SFM_DATA_HPP_
#define INCLUDE_EDGEGRAPH3D_IO_OUTPUT_OUTPUT_SFM_DATA_HPP_

#include <string>
#include <vector>

#include "global_defines.hpp"
#include "global_switches.hpp"

struct SfMData;

void output_sfm_data(const char *original_sfmd_filename, const SfMData &sfmd, const std::string &output_file);
void output_filtered_sfm_data(const char *original_sfmd_filename, const SfMData &sfmd, const std::vector<bool> &inliers, const std::string &output_file);

#endif /* INCLUDE_EDGEGRAPH3D_IO_OUTPUT_OUTPUT_SFM_DATA_HPP_ */
