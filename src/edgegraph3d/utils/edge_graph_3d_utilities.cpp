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


#include "edge_graph_3d_utilities.hpp"

#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <opencv2/core/base.hpp>
#include <opencv2/core/cvdef.h>
#include <opencv2/core/cvstd.inl.hpp>
#include <opencv2/core/fast_math.hpp>
#include <opencv2/core/hal/interface.h>
#include <opencv2/core/mat.inl.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs/imgcodecs_c.h>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>

#include <stddef.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>

#include <glm.hpp>
#include "OpenMvgParser.h"
#include "SfMData.h"

#include "geometric_utilities.hpp"

using namespace std;
using namespace cv;

void print_vec2(const glm::vec2 &p2d) {
	cout << "[" << p2d[0] << "," << p2d[1] << "]";
}

void print_vec4(const glm::vec4 &p4d) {
	cout << "[" << p4d[0] << "," << p4d[1] << "," << p4d[2] << "," << p4d[3] << "]";
}

void print_vector_vec2(const vector<glm::vec2> &p2ds) {
	cout << "[";
	if(p2ds.size() > 0) {
		print_vec2(p2ds[0]);
		for(int i=1;i<p2ds.size();i++) {
			cout << ",";
			print_vec2(p2ds[i]);
		}
	}
	cout << "]";
}

float avg_vec_float(const vector<float> &v) {
	if(v.empty())
		return 0;

	float avg=0;
	for(const auto el: v)
		avg += el;
	avg /= v.size();

	return avg;
}

float avg_vec_float_weight_starting_image(const vector<float> &v, const vector<int> &selected_id, const int starting_img, const float starting_weight) {
	if(v.empty())
		return 0;

	float avg=0;
	for(int i=0; i < v.size(); i++) {
		if(selected_id[i] != starting_img)
			avg += v[i];
		else
			avg += (v[i]*starting_weight);
	}
	avg /= v.size();

	return avg;
}

vector<Mat> copy_vector_Mat(vector<Mat> &matvec) {
	vector<Mat> vm;
	for(auto m : matvec)
		vm.push_back(m.clone());
	return vm;
}

void release_vector_Mat(vector<Mat> &matvec) {
	for(auto &m : matvec)
		m.release();
	matvec.clear();
}

void print_glmmat3(glm::mat3 m) {
	const int r = 3;
	const int c = 3;

	cout << "[";
	for(int i=0; i < r; i++) {
		cout << "\t" << m[i][0];
		for(int j=1; j < c; j++)
			cout << ",\t" << m[i][j];
		if(i < r-1)
			cout << endl;
	}
	cout << "]" << endl;
}

void print_glmmat4(glm::mat4 m) {
	const int r = 4;
	const int c = 4;

	cout << "[";
	for(int i=0; i < r; i++) {
		cout << "\t" << m[i][0];
		for(int j=1; j < c; j++)
			cout << ",\t" << m[i][j];
		if(i < r-1)
			cout << endl;
	}
	cout << "]" << endl;
}

float rand_float() {
	return static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
}

vector<Mat> copy_imgs(const vector<Mat> &imgs) {
	vector<Mat> res;
	for(const auto &i: imgs)
		res.push_back(i.clone());
	return res;
}

vector<glm::vec4> convert_Vec4f_vec4(const vector<Vec4f> &cur_edges) {
	vector<glm::vec4> glm_edges;

	for (auto& v: cur_edges)
		glm_edges.push_back(glm::vec4(v[0],v[1],v[2],v[3]));

	return glm_edges;
}

vector<Vec4f> convert_vec4_Vec4f(const vector<glm::vec4> &cur_edges) {
	vector<Vec4f> edges_Vec4f;

	for (auto& e: cur_edges) {
		Vec4f v(e[0],e[1],e[2],e[3]);

		edges_Vec4f.push_back(v);
	}

	return edges_Vec4f;
}



void convert_glm_mat4_to_cv_Mat34(const glm::mat4 &glm_mat, cv::Mat &out) {
	out = cv::Mat(3, 4, CV_32F);
	for (int row = 0; row < 3; row++) {
		for (int col = 0; col < 4; col++) {
			out.at<float>(row, col) = glm_mat[row][col];
		}
	}
}

void convert_glm_mat4_to_cv_Mat(const glm::mat4 &glm_mat, cv::Mat &out) {
	out = cv::Mat(4, 4, CV_32F);
	for (int row = 0; row < 4; row++) {
		for (int col = 0; col < 4; col++) {
			out.at<float>(row, col) = glm_mat[row][col];
		}
	}
}

// http://stackoverflow.com/questions/13299409/how-to-get-the-image-pixel-at-real-locations-in-opencv
cv::Vec3b getColorSubpix(const cv::Mat& img, cv::Point2f pt)
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

// Returns color of given point, by averaging color of all its 2D observations
Scalar get_refpoint_mixed_color(const SfMData &sfmd, const vector<Mat> &imgs, const ulong refpoint_id) {
	Scalar out_col(0,0,0);
	int numviews = sfmd.camViewingPointN_[refpoint_id].size();
	for(int i=0; i < numviews; i++)
		out_col += Scalar(getColorSubpix(imgs[sfmd.camViewingPointN_[refpoint_id][i]],convert_glm_vec2_to_cv_Point2f(sfmd.point2DoncamViewingPoint_[refpoint_id][i])));
	return Scalar(out_col[0]/float(numviews),out_col[1]/float(numviews),out_col[2]/float(numviews));
}

cv::Point2f convert_glm_vec2_to_cv_Point2f(const glm::vec2 &glm_vec) {
	cv::Point2f out;
	out.x = glm_vec[0];
	out.y = glm_vec[1];
	return out;
}

cv::Point3f convert_glm_vec3_to_cv_Point3f(const glm::vec3 &glm_vec) {
	return cv::Point3f(glm_vec[0],glm_vec[1],glm_vec[2]);
}

void convert_glm_vec2_to_cv_Point2f(const glm::vec2 &glm_vec, cv::Point2f &out) {
	out.x = glm_vec[0];
	out.y = glm_vec[1];
}

void convert_glm_vec2_to_cv_Mat21(const glm::vec2 &glm_vec, cv::Mat &out) {
	out = cv::Mat(2, 1, CV_64F);
	out.at<double>(0, 0) = glm_vec[0];
	out.at<double>(1, 0) = glm_vec[1];
}

std::string remove_path(std::string s) {
	std::size_t found = s.find_last_of("/\\");
	return s.substr(found+1);
}

std::string remove_extension(std::string s) {
	std::size_t found = s.find_last_of(".");
	return s.substr(0,found);
}


std::string remove_path_and_exception(std::string s) {
	std::size_t found = s.find_last_of("/\\");
	size_t lastindex = s.find_last_of(".");
	return s.substr(found+1,lastindex);
}



vector<Mat> parse_images(const char *images_folder, const SfMData &sfm_data) {
	vector<Mat> imgs;
	Mat img;
	string img_path;

	for(auto s : sfm_data.camerasPaths_) {
		  img_path = images_folder + remove_path(s);
		  img = imread(img_path, CV_LOAD_IMAGE_COLOR);

		  if(!img.data) {
			  cout << ("Could not read \"" + img_path + "\"\n");
			  std::invalid_argument("Input image not found");
		  }

		  imgs.push_back(img);
	}

	return imgs;
}

vector<Mat> parse_imagesRGB(const char *images_folder, const SfMData &sfm_data) {
	vector<Mat> imgs;
	Mat img;
	string img_path;

	for(auto s : sfm_data.camerasPaths_) {
		  img_path = images_folder + remove_path(s);
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

int parse_images(const char *images_folder, const SfMData &sfm_data, vector<Mat> &imgs) {
	Mat img;
	string img_path;
	imgs.clear();

	for(auto s : sfm_data.camerasPaths_) {
		  img_path = images_folder + remove_path(s);
		  img = imread(img_path, CV_LOAD_IMAGE_COLOR);

		  if(!img.data) {
			  cout << ("Could not read \"" + img_path + "\"\n");
			  return -1;
		  }

		  imgs.push_back(img);
	}

	return 0;
}

set<int> find_points_on_both_images(const vector<set<int>> &points_on_images, int img_i, int img_j) {
	const set<int> &s_i = points_on_images[img_i];
	const set<int> &s_j = points_on_images[img_j];
	set<int> points_on_both;
	set_intersection(s_i.begin(),s_i.end(),s_j.begin(),s_j.end(),
		                  std::inserter(points_on_both,points_on_both.begin()));
	return points_on_both;
}

vector<pair<glm::vec2,glm::vec2>> find_fakepoints_on_both_images(const SfMData &sfmd, int img_i, int img_j, const ulong amount) {
	vector<pair<glm::vec2,glm::vec2>> res;

	ulong actual_amount = sfmd.pointsVisibleFromCamN_[img_i].size() > amount ? amount : sfmd.pointsVisibleFromCamN_[img_i].size();

	for(ulong i=0; i < actual_amount; i++) {
		ulong refpoint_id = sfmd.pointsVisibleFromCamN_[img_i][i];
		glm::vec2 refpoint_coords1 = sfmd.point2DoncamViewingPoint_[img_i][i];
		glm::vec2 refpoint_coords2 = compute_projection(sfmd,img_j,sfmd.points_[refpoint_id]);
		res.push_back(make_pair(refpoint_coords1,refpoint_coords2));
	}

	return res;
}

vector<set<int>> get_point_sets_on_images(const SfMData &sfmd) {
	vector<set<int>> points_on_images;
	set<int> cur_points;

	for(auto &points_on_img: sfmd.pointsVisibleFromCamN_) {
		cur_points.clear();
		std::copy( points_on_img.begin(), points_on_img.end(), std::inserter( cur_points, cur_points.end() ) );
		points_on_images.push_back(cur_points);
	}

	return points_on_images;
}

glm::vec2 get_2d_coordinates_of_point_on_image(const SfMData &sfmd, const int img_id, const int point_id, bool &found) {
	found = false;
	glm::vec2 point;

	for(int i=0; i < sfmd.camViewingPointN_[point_id].size(); i++)
		if(sfmd.camViewingPointN_[point_id][i] == img_id) {
			point = sfmd.point2DoncamViewingPoint_[point_id][i];
			found = true;
		}

	return point;
}

void get_2d_coordinates_of_point_on_image(const SfMData &sfmd, const int img_id, const int point_id, glm::vec2 &point, bool &found) {
	found = false;

	for(int i=0; i < sfmd.camViewingPointN_[point_id].size(); i++)
		if(sfmd.camViewingPointN_[point_id][i] == img_id) {
			point = sfmd.point2DoncamViewingPoint_[point_id][i];
			found = true;
		}
}

vector<glm::vec2> get_vector_of_2d_coordinates_of_point_on_image(const SfMData &sfmd, const int img_id, const vector<int> &point_ids) {
	bool found;
	glm::vec2 p;
	vector<glm::vec2> result;

	for(const auto point_id : point_ids) {
		get_2d_coordinates_of_point_on_image(sfmd, img_id, point_id, p, found);
		result.push_back(p);
	}

	return result;
}

cv::Vec3f convert_from_glm_vec3_to_cv_Vec3f(const glm::vec3 &v, cv::Vec3f &vcv) {
	vcv[0]=v[0];
	vcv[1]=v[1];
	vcv[2]=v[2];
}

cv::Vec3f convert_from_glm_vec3_to_cv_Mat(const glm::vec3 &v, cv::Mat &vcv) {
	vcv.at<float>(0,0)=v[0];
	vcv.at<float>(1,0)=v[1];
	vcv.at<float>(2,0)=v[2];
}

void convert_from_glm_mat3_to_cv_Mat3f(const glm::mat3 &m, cv::Mat &mcv) {
	mcv.at<float>(0,0)=m[0][0];
	mcv.at<float>(0,1)=m[0][1];
	mcv.at<float>(0,2)=m[0][2];
	mcv.at<float>(1,0)=m[1][0];
	mcv.at<float>(1,1)=m[1][1];
	mcv.at<float>(1,2)=m[1][2];
	mcv.at<float>(2,0)=m[2][0];
	mcv.at<float>(2,1)=m[2][1];
	mcv.at<float>(2,2)=m[2][2];
}

int find_closest_point_on_image(const SfMData &sfmd, const int img_id, const glm::vec2 &p) {
	float min_dist,cur_dist;
	int cur_min;
	bool found;
	glm::vec2 other_p;
	const std::vector<int> &visible_points = sfmd.pointsVisibleFromCamN_[img_id];

	get_2d_coordinates_of_point_on_image(sfmd,img_id,visible_points[0],other_p,found);
	min_dist = squared_2d_distance(p,other_p);
	cur_min = visible_points[0];

	for(const auto other_p_id : visible_points) {
		get_2d_coordinates_of_point_on_image(sfmd,img_id,other_p_id,other_p,found);
		cur_dist = squared_2d_distance(p,other_p);
		// cout << "Point " << other_p_id << " -> dist : " << cur_dist << endl;
		if(cur_dist < min_dist) {
			cur_min = other_p_id;
			min_dist = cur_dist;
		}
	}

	return cur_min;
}

bool glm_vec2_equal(const glm::vec2 &v1, const glm::vec2 &v2) {
	return v1 == v2;
}

bool vec_glm_vec2_equal(const vector<glm::vec2> &v1, const vector<glm::vec2> &v2) {
	if(v1.size() != v2.size())
		return false;

	for(ulong i = 0; i < v1.size(); i++)
		if(!glm_vec2_equal(v1[i],v2[i]))
			return false;

	return true;
}

bool glm_vec3_equal(const glm::vec3 &v1, const glm::vec3 &v2) {
	return v1 == v2;
}

bool vec_glm_vec3_equal(const vector<glm::vec3> &v1, const vector<glm::vec3> &v2) {
	if(v1.size() != v2.size())
		return false;

	for(ulong i = 0; i < v1.size(); i++)
		if(!glm_vec3_equal(v1[i],v2[i]))
			return false;

	return true;
}

// check if two vector of vec2 are equal but reversed
bool vec_glm_vec2_equal_inv(const vector<glm::vec2> &v1, const vector<glm::vec2> &v2) {
	if(v1.size() != v2.size())
		return false;

	for(ulong i = 0; i < v1.size(); i++)
		if(!glm_vec2_equal(v1[i],v2[v1.size()-i-1]))
			return false;

	return true;
}

bool vec_glm_vec3_equal_inv(const vector<glm::vec3> &v1, const vector<glm::vec3> &v2) {
	if(v1.size() != v2.size())
		return false;

	for(ulong i = 0; i < v1.size(); i++)
		if(!glm_vec3_equal(v1[i],v2[v1.size()-i-1]))
			return false;

	return true;
}

std::ostream &operator<< (std::ostream &out, const glm::vec2 &vec) {
    out << "("
        << vec.x << ", " << vec.y
        << ")";

    return out;
}

std::ostream &operator<< (std::ostream &out, const glm::vec3 &vec) {
    out << "("
        << vec.x << ", " << vec.y << ", "<< vec.z
        << ")";

    return out;
}

std::ostream &operator<< (std::ostream &out, const glm::vec4 &vec) {
    out << "("
        << vec.x << ", " << vec.y << ", "<< vec.z << ", "<< vec.w
        << ")";

    return out;
}

std::vector<glm::vec2> intersections_remove_segments(const std::vector<SegmentEdgeManager::segment_point> &intersections) {
	std::vector<glm::vec2> res;
	for(const auto &i : intersections)
		res.push_back(i.coords);
	return res;
}

std::vector<std::vector<std::vector<glm::vec2>>> correspondences_remove_segments(const std::vector<std::vector<std::vector<SegmentEdgeManager::segment_point>>> &correspondences) {
	std::vector<std::vector<std::vector<glm::vec2>>> res;
	for(const auto &v1 : correspondences) {
		std::vector<std::vector<glm::vec2>> v1c;
		for(const auto &v2 : v1) {
			std::vector<glm::vec2> v2c;
			for(const auto &corresp : v2)
				v2c.push_back(corresp.coords);
			v1c.push_back(v2c);
		}
		res.push_back(v1c);
	}
	return res;
}

glm::mat4x4 get_rt4x4(const glm::mat3 &r, const glm::vec3 &t) {
	return glm::mat4(
			r[0][0],r[0][1],r[0][2],t[0],
			r[1][0],r[1][1],r[1][2],t[1],
			r[2][0],r[2][1],r[2][2],t[2],
			0.0,	0.0,	0.0,	1.0);
}

std::ostream &operator<<(std::ostream &out, const glm::mat3 &m) {
	return out << "[" << m[0][0] << "," << m[0][1] << "," << m[0][2] << ";" << m[1][0] << "," << m[1][1] << "," << m[1][2] << ";" << m[2][0] << "," << m[2][1] << "," << m[2][2] << "]";
}

std::ostream &operator<<(std::ostream &out, const glm::mat4 &m) {
	return out << "[" << m[0][0] << "," << m[0][1] << "," << m[0][2] << ","  << m[0][3] << ";" << m[1][0] << "," << m[1][1] << "," << m[1][2] << ","  << m[1][3] << ";" << m[2][0] << "," << m[2][1] << "," << m[2][2] << ","  << m[2][3] << ";" << m[3][0] << "," << m[3][1] << "," << m[3][2] << "," << m[3][3] << "]";
}

Mat** create_2D_Mat_array(unsigned long nrows, unsigned long ncols) {
	Mat** ptr = create_2D_array<Mat>(nrows, ncols);

	for(int i=0; i < nrows;i++)
		for(int j=0; j < nrows;j++)
			ptr[i][j] = Mat::zeros(3, 3, CV_32FC2);

	return ptr;
}

void delete_2D_Mat_array(Mat** arr, unsigned long nrows, unsigned long ncols) {
	for(int i=0; i < nrows;i++)
		for(int j=0; j < nrows;j++) {
			arr[i][j].release();
		}

	delete_2D_array<Mat>(arr);
}

#define SMALL_FLOAT_EPSILON 0.001

float floor_or_upper_if_close(const float v) {
	if(ceil(v) - v < SMALL_FLOAT_EPSILON)
		return ceil(v);
	else
		return floor(v);
}

bool is_m_multiple_of_n_float(const float m, const float n) {
	float div = m / n;
	float mul = floor_or_upper_if_close(div) * n;
	if(abs(m - mul) < SMALL_FLOAT_EPSILON)
		return true;
	else
		return false;
}

pair<ulong,ulong> get_2dmap_cell_from_coords(const float cell_dim, const glm::vec2 &coords, bool &on_boundary) {
	on_boundary = is_m_multiple_of_n_float(coords.x,cell_dim) || is_m_multiple_of_n_float(coords.y,cell_dim);

	return make_pair((int) floor_or_upper_if_close(coords.x / cell_dim), (int) floor_or_upper_if_close(coords.y / cell_dim));
}

pair<ulong,ulong> get_2dmap_cell_from_coords(const float cell_dim, const glm::vec2 &coords, bool &on_boundary_row, bool &on_boundary_col)  {
	on_boundary_row = is_m_multiple_of_n_float(coords.x,cell_dim);
	on_boundary_col = is_m_multiple_of_n_float(coords.y,cell_dim);

	return make_pair((int) floor_or_upper_if_close(coords.x / cell_dim), (int) floor_or_upper_if_close(coords.y / cell_dim));
}

string compute_plg_path(char *out_images_folder, const SfMData &sfm_data, int img_id, const String &s) {
	return out_images_folder + s + remove_path_and_exception(sfm_data.camerasPaths_[img_id]) + PLG_EXTENSION;
}

string compute_img_path(char *out_images_folder, const SfMData &sfm_data, int img_id, const String &s) {
	return out_images_folder + s + remove_path(sfm_data.camerasPaths_[img_id]);
}

void write_img(const string &out_images_folder, const SfMData &sfm_data, int img_id, const Mat &img, const String &s) {
	string img_path = compute_img_path(out_images_folder.c_str(),sfm_data,img_id,s);
	cout << "Writing to " << img_path << endl;
	imwrite(img_path,img);
}

int write_images(const string &out_images_folder, const SfMData &sfm_data, vector<Mat> &imgs, const String &s) {

	cout << "Saving images:" << endl;

	for(int i=0; i < imgs.size(); i++)
		write_img(out_images_folder,sfm_data,i,imgs[i],s);

	return 0;
}

int write_images(const string &out_images_folder, const SfMData &sfm_data, vector<Mat> &imgs) {
	return write_images(out_images_folder, sfm_data, imgs, "");
}

pair<vector<glm::vec3>,vector<Scalar>> parseGTxyz(const string &spath) {
	const boost::filesystem::path path(spath);
	vector<glm::vec3> gtpointcloud_;
	vector<Scalar> gtpointcloudcolors_;
	double x,y,z,cx,cy,cz;

	if(!(boost::filesystem::exists(path) && boost::filesystem::is_regular_file(path) && boost::filesystem::extension(path) == ".xyz"))
		std::invalid_argument("Cannot parge groundtruth cloud");

	std::ifstream f(path.c_str());

	if(!f.is_open())
		std::invalid_argument("Cannot parge groundtruth cloud");

	std::string line;
	float alpha;
	while (std::getline(f, line)) {
		std::istringstream ss(line);
		ss >> x >> y >> z >> cx >> cy >> cz >> alpha;
		gtpointcloud_.push_back(glm::vec3(x,y,z));
		gtpointcloudcolors_.push_back(Scalar(cz,cy,cx));
	}

	return make_pair(gtpointcloud_, gtpointcloudcolors_);
}

pair<vector<glm::vec3>,vector<Scalar>> parseGTply(const string &spath) {
	const boost::filesystem::path path(spath);
	vector<glm::vec3> gtpointcloud_;
	vector<Scalar> gtpointcloudcolors_;
	double x,y,z,cx,cy,cz;

	if(!(boost::filesystem::exists(path) && boost::filesystem::is_regular_file(path) && boost::filesystem::extension(path) == ".xyz"))
		std::invalid_argument("Cannot parge groundtruth cloud");

	std::ifstream f(path.c_str());

	if(!f.is_open())
		std::invalid_argument("Cannot parge groundtruth cloud");

	std::string line;
	float alpha;

	bool found_ply_header_end=false;

	while (std::getline(f, line))
		if(strcmp(line.c_str(),"end_header") == 0) {
			found_ply_header_end=true;
			break;
		}

	if(!found_ply_header_end)
		std::invalid_argument("Cannot parge groundtruth cloud - PLY header not found");

	while (std::getline(f, line)) {
		std::istringstream ss(line);
		ss >> x >> y >> z >> cx >> cy >> cz >> alpha;
		gtpointcloud_.push_back(glm::vec3(x,y,z));
		gtpointcloudcolors_.push_back(Scalar(cz,cy,cx));
	}

	return make_pair(gtpointcloud_, gtpointcloudcolors_);
}

glm::vec3 to_glmvec3(CGAL::Simple_cartesian<double>::Point_3 p) {
	return glm::vec3(p[0],p[1],p[2]);
}

std::ostream &operator<< (std::ostream &out, const pair<vector<int>,vector<ulong>> &p)
{
	out << "{" << p.first << "," << p.second << "}";
	return out;
}

SfMData read_sfm_data(const char *sfmd_path) {
	  OpenMvgParser op(sfmd_path);
	  op.parse();
	  return op.getSfmData();
}








