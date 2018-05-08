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


#ifndef INCLUDE_UTILS_OUTPUT_POINT_CLOUD_HPP_
#define INCLUDE_UTILS_OUTPUT_POINT_CLOUD_HPP_

#include <opencv2/core/types.hpp>
#include <string>
#include <vector>

#include "global_defines.hpp"
#include "global_switches.hpp"

struct SfMData;

void output_point_cloud(const SfMData &sfmd, const char *out_folder, const std::string &output_filename);
void output_point_cloud(const SfMData &sfmd, const char *output_filename);

void output_colored_point_cloud(const SfMData &sfmd, const char *images_folder, const char *out_folder, const std::string &output_filename);

void output_point_cloud(const SfMData &sfmd, const std::vector<bool> &inliers, const char *out_folder, const std::string &output_filename);

void output_colored_point_cloud(const SfMData &sfmd, const char *images_folder, const char *out_folder, const std::string &output_filename);

void output_colored_point_cloud(const SfMData &sfmd, const char *images_folder, const char *output_filename);

// Draws inliers green, outliers red
void output_colored_point_cloud_inliers_highlight(const SfMData &sfmd, const std::vector<bool> &inliers, const cv::Scalar &inliers_color, const cv::Scalar &outliers_color, const char *out_folder, const std::string &output_filename);
void output_colored_point_cloud_2_inliers_vec_highlight(const SfMData &sfmd, const std::vector<bool> &old_inliers, const std::vector<bool> &new_inliers, const cv::Scalar &inliers_color, const cv::Scalar &outliers_color, const cv::Scalar &new_outliers_color, const char *out_folder, const std::string &output_filename);
void output_colored_point_cloud_sor_inliers_highlight(const SfMData &sfmd, const std::vector<bool> &inliers, const std::vector<bool> &sor_inliers, const cv::Scalar &inliers_color, const cv::Scalar &outliers_color, const cv::Scalar &sor_outliers_color, const char *out_folder, const std::string &output_filename);

#endif /* INCLUDE_UTILS_OUTPUT_POINT_CLOUD_HPP_ */
