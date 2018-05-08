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


#include "output_point_cloud.hpp"

#include <opencv2/core/base.hpp>
#include <opencv2/core/fast_math.hpp>
#include <opencv2/core/hal/interface.h>
#include <opencv2/core/mat.hpp>
#include <opencv2/core/mat.inl.hpp>
#include <opencv2/core/matx.hpp>
#include <opencv2/core/types.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs/imgcodecs_c.h>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <sys/types.h>
#include <cassert>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "SfMData.h"
#include "glm.hpp"
#include "glm.hpp"

using namespace std;
using namespace cv;

std::string opc_remove_path(std::string s) {
	std::size_t found = s.find_last_of("/\\");
	return s.substr(found+1);
}

cv::Vec3b opc_getColorSubpix(const cv::Mat& img, cv::Point2f pt)
{
    assert(!img.empty());
    assert(img.channels() == 3);

    int x = (int)pt.x;
    int y = (int)pt.y;

    int x0 = cv::borderInterpolate(x,   img.cols, cv::BORDER_REFLECT_101);
    int x1 = cv::borderInterpolate(x+1, img.cols, cv::BORDER_REFLECT_101);
    int y0 = cv::borderInterpolate(y,   img.rows, cv::BORDER_REFLECT_101);
    int y1 = cv::borderInterpolate(y+1, img.rows, cv::BORDER_REFLECT_101);

    float a = pt.x - (float)x;
    float c = pt.y - (float)y;

    uchar b = (uchar)cvRound((img.at<cv::Vec3b>(y0, x0)[0] * (1.f - a) + img.at<cv::Vec3b>(y0, x1)[0] * a) * (1.f - c)
                           + (img.at<cv::Vec3b>(y1, x0)[0] * (1.f - a) + img.at<cv::Vec3b>(y1, x1)[0] * a) * c);
    uchar g = (uchar)cvRound((img.at<cv::Vec3b>(y0, x0)[1] * (1.f - a) + img.at<cv::Vec3b>(y0, x1)[1] * a) * (1.f - c)
                           + (img.at<cv::Vec3b>(y1, x0)[1] * (1.f - a) + img.at<cv::Vec3b>(y1, x1)[1] * a) * c);
    uchar r = (uchar)cvRound((img.at<cv::Vec3b>(y0, x0)[2] * (1.f - a) + img.at<cv::Vec3b>(y0, x1)[2] * a) * (1.f - c)
                           + (img.at<cv::Vec3b>(y1, x0)[2] * (1.f - a) + img.at<cv::Vec3b>(y1, x1)[2] * a) * c);

    return cv::Vec3b(b, g, r);
}

cv::Point2f opc_convert_glm_vec2_to_cv_Point2f(const glm::vec2 &glm_vec) {
	cv::Point2f out;
	out.x = glm_vec[0];
	out.y = glm_vec[1];
	return out;
}

// Returns color of given point, by averaging color of all its 2D observations
Scalar opc_get_refpoint_mixed_color(const SfMData &sfmd, const vector<Mat> &imgs, const ulong refpoint_id) {
	Scalar out_col(0,0,0);
	int numviews = sfmd.camViewingPointN_[refpoint_id].size();
	for(int i=0; i < numviews; i++)
		out_col += Scalar(opc_getColorSubpix(imgs[sfmd.camViewingPointN_[refpoint_id][i]],opc_convert_glm_vec2_to_cv_Point2f(sfmd.point2DoncamViewingPoint_[refpoint_id][i])));
	return Scalar(out_col[0]/float(numviews),out_col[1]/float(numviews),out_col[2]/float(numviews));
}

vector<Mat> opc_parse_imagesRGB(const char *images_folder, const SfMData &sfm_data) {
	vector<Mat> imgs;
	Mat img;
	string img_path;

	for(auto s : sfm_data.camerasPaths_) {
		  img_path = images_folder + opc_remove_path(s);
		  img = imread(img_path, CV_LOAD_IMAGE_COLOR);

		  if(!img.data) {
			  cout << ("Could not read \"" + img_path + "\"\n");
			  std::invalid_argument("Input image not found");
		  }
		  cv::cvtColor(img, img, COLOR_BGR2RGB);
		  imgs.push_back(img);
	}

	return imgs;
}

void output_point_cloud(const SfMData &sfmd, const char *out_folder, const string &output_filename) {
	output_point_cloud(sfmd, (std::string(out_folder)+output_filename).c_str());
}

void output_point_cloud(const SfMData &sfmd, const char *output_filename) {
	ofstream outfile;
	outfile.open(output_filename);

	const Scalar whitecol(255,255,255);

	outfile <<
"ply" << endl <<
"format ascii 1.0" << endl <<
"element vertex " << sfmd.numPoints_ << endl <<
"property float x" << endl <<
"property float y" << endl <<
"property float z" << endl <<
"property uchar red" << endl <<
"property uchar green" << endl <<
"property uchar blue" << endl <<
"end_header" << endl;

	for(ulong i=0; i < sfmd.numPoints_; i++) {
		outfile << sfmd.points_[i].x << " " << sfmd.points_[i].y << " " << sfmd.points_[i].z << " " << int(whitecol[0]) << " " << int(whitecol[1]) << " " << int(whitecol[2]);
		if(i<sfmd.numPoints_-1)
			outfile << endl;
	}

	outfile.close();
}

void output_point_cloud(const SfMData &sfmd, const vector<bool> &inliers, const char *out_folder, const string &output_filename) {
	ofstream outfile;
	outfile.open((std::string(out_folder)+output_filename).c_str());

	const Scalar whitecol(255,255,255);

	ulong npoints=0;
	for(ulong i=0; i < inliers.size();i++)
		if(inliers[i])
			npoints++;

	outfile <<
"ply" << endl <<
"format ascii 1.0" << endl <<
"element vertex " << npoints << endl <<
"property float x" << endl <<
"property float y" << endl <<
"property float z" << endl <<
"property uchar red" << endl <<
"property uchar green" << endl <<
"property uchar blue" << endl <<
"end_header" << endl;

	for(ulong i=0; i < sfmd.numPoints_; i++)
		if(inliers[i]) {
			outfile << sfmd.points_[i].x << " " << sfmd.points_[i].y << " " << sfmd.points_[i].z << " " << int(whitecol[0]) << " " << int(whitecol[1]) << " " << int(whitecol[2]);
			if(i<sfmd.numPoints_-1)
				outfile << endl;
		}

	outfile.close();
}

void output_colored_point_cloud(const SfMData &sfmd, const char *images_folder, const char *out_folder, const string &output_filename) {
	output_point_cloud(sfmd, images_folder, (std::string(out_folder)+output_filename).c_str());
}

void output_colored_point_cloud(const SfMData &sfmd, const char *images_folder, const char *output_filename) {
	vector<Mat> imgs = opc_parse_imagesRGB(images_folder,sfmd);
	ofstream outfile;
	outfile.open(output_filename);

	Scalar cur_color;

	outfile <<
"ply" << endl <<
"format ascii 1.0" << endl <<
"element vertex " << sfmd.numPoints_ << endl <<
"property float x" << endl <<
"property float y" << endl <<
"property float z" << endl <<
"property uchar red" << endl <<
"property uchar green" << endl <<
"property uchar blue" << endl <<
"end_header" << endl;

	for(ulong i=0; i < sfmd.numPoints_; i++) {
		cur_color = opc_get_refpoint_mixed_color(sfmd, imgs, i);
		outfile << sfmd.points_[i].x << " " << sfmd.points_[i].y << " " << sfmd.points_[i].z << " " << int(cur_color[0]) << " " << int(cur_color[1]) << " " << int(cur_color[2]);
		if(i<sfmd.numPoints_-1)
			outfile << endl;
	}

	outfile.close();
}

void output_colored_point_cloud(const SfMData &sfmd, const char *images_folder, const vector<bool> &inliers, const char *out_folder, const string &output_filename) {
	vector<Mat> imgs = opc_parse_imagesRGB(images_folder,sfmd);
	ofstream outfile;
	outfile.open(out_folder+output_filename);

	Scalar cur_color;

	ulong npoints=0;
	for(ulong i=0; i < inliers.size();i++)
		if(inliers[i])
			npoints++;

	outfile <<
"ply" << endl <<
"format ascii 1.0" << endl <<
"element vertex " << npoints << endl <<
"property float x" << endl <<
"property float y" << endl <<
"property float z" << endl <<
"property uchar red" << endl <<
"property uchar green" << endl <<
"property uchar blue" << endl <<
"end_header" << endl;

	for(ulong i=0; i < sfmd.numPoints_; i++)
		if(inliers[i]) {
			cur_color = opc_get_refpoint_mixed_color(sfmd, imgs, i);
			outfile << sfmd.points_[i].x << " " << sfmd.points_[i].y << " " << sfmd.points_[i].z << " " << int(cur_color[0]) << " " << int(cur_color[1]) << " " << int(cur_color[2]);
			if(i<sfmd.numPoints_-1)
				outfile << endl;
		}

	outfile.close();
}

// Draws inliers green, outliers red
void output_colored_point_cloud_inliers_highlight(const SfMData &sfmd, const vector<bool> &inliers, const Scalar &inliers_color, const Scalar &outliers_color, const char *out_folder, const string &output_filename) {
	ofstream outfile;
	outfile.open(out_folder+output_filename);

	outfile <<
"ply" << endl <<
"format ascii 1.0" << endl <<
"element vertex " << sfmd.numPoints_ << endl <<
"property float x" << endl <<
"property float y" << endl <<
"property float z" << endl <<
"property uchar red" << endl <<
"property uchar green" << endl <<
"property uchar blue" << endl <<
"end_header" << endl;

	for(ulong i=0; i < sfmd.numPoints_; i++){
			const Scalar &cur_color = inliers[i] ? inliers_color : outliers_color;
			outfile << sfmd.points_[i].x << " " << sfmd.points_[i].y << " " << sfmd.points_[i].z << " " << int(cur_color[0]) << " " << int(cur_color[1]) << " " << int(cur_color[2]);
			if(i<sfmd.numPoints_-1)
				outfile << endl;
		}

	outfile.close();
}

// Draws inliers green, outliers red, outliers_sor orange
void output_colored_point_cloud_2_inliers_vec_highlight(const SfMData &sfmd, const vector<bool> &old_inliers, const vector<bool> &new_inliers, const Scalar &inliers_color, const Scalar &outliers_color, const Scalar &new_outliers_color, const char *out_folder, const string &output_filename) {
	ofstream outfile;
	outfile.open(out_folder+output_filename);

	outfile <<
"ply" << endl <<
"format ascii 1.0" << endl <<
"element vertex " << sfmd.numPoints_ << endl <<
"property float x" << endl <<
"property float y" << endl <<
"property float z" << endl <<
"property uchar red" << endl <<
"property uchar green" << endl <<
"property uchar blue" << endl <<
"end_header" << endl;

	//Scalar cur_color;
	for(ulong i=0; i < sfmd.numPoints_; i++){
			const Scalar &cur_color = !old_inliers[i] ? outliers_color : (new_inliers[i] ? inliers_color : new_outliers_color);

			outfile << sfmd.points_[i].x << " " << sfmd.points_[i].y << " " << sfmd.points_[i].z << " " << int(cur_color[0]) << " " << int(cur_color[1]) << " " << int(cur_color[2]);
			if(i<sfmd.numPoints_-1)
				outfile << endl;
		}

	outfile.close();
}

// Draws inliers green, outliers red, outliers_sor orange
void output_colored_point_cloud_sor_inliers_highlight(const SfMData &sfmd, const vector<bool> &inliers, const vector<bool> &sor_inliers, const Scalar &inliers_color, const Scalar &outliers_color, const Scalar &sor_outliers_color, const char *out_folder, const string &output_filename) {
	ofstream outfile;
	outfile.open(out_folder+output_filename);

	outfile <<
"ply" << endl <<
"format ascii 1.0" << endl <<
"element vertex " << sfmd.numPoints_ << endl <<
"property float x" << endl <<
"property float y" << endl <<
"property float z" << endl <<
"property uchar red" << endl <<
"property uchar green" << endl <<
"property uchar blue" << endl <<
"end_header" << endl;

	int sor_index=0;
	//Scalar cur_color;
	for(ulong i=0; i < sfmd.numPoints_; i++){
			const Scalar &cur_color = !inliers[i] ? outliers_color : (sor_inliers[sor_index] ? inliers_color : sor_outliers_color);

			if(inliers[i])
				sor_index++;

			outfile << sfmd.points_[i].x << " " << sfmd.points_[i].y << " " << sfmd.points_[i].z << " " << int(cur_color[0]) << " " << int(cur_color[1]) << " " << int(cur_color[2]);
			if(i<sfmd.numPoints_-1)
				outfile << endl;
		}

	outfile.close();
}


