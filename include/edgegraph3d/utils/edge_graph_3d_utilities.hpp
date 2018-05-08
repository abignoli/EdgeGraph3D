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


#ifndef INCLUDE_EDGEGRAPH3D_UTILS_EDGE_GRAPH_3D_UTILITIES_HPP_
#define INCLUDE_EDGEGRAPH3D_UTILS_EDGE_GRAPH_3D_UTILITIES_HPP_

#include <CGAL/Simple_cartesian.h>
#include <opencv2/core/cvstd.hpp>
#include <opencv2/core/mat.hpp>
#include <opencv2/core/matx.hpp>
#include <opencv2/core/types.hpp>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "segment_edge_manager.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"
#include "glm.hpp"


#define PLG_EXTENSION ".plg"

#define DO_MEASURE_OMP_TIME_START double DO_MEASURE_OMP_TIME_START_start = omp_get_wtime();
#define DO_MEASURE_OMP_TIME_END double DO_MEASURE_OMP_TIME_END_end = omp_get_wtime(); cout << "DO_MEASURE_OMP_TIME: " << DO_MEASURE_OMP_TIME_END_end - DO_MEASURE_OMP_TIME_START_start << "\n";
#define DO_MEASURE_OMP_TIME_END_NAMED(X) double DO_MEASURE_OMP_TIME_END_end = omp_get_wtime(); cout << "DO_MEASURE_OMP_TIME (" << X << "): " << DO_MEASURE_OMP_TIME_END_end - DO_MEASURE_OMP_TIME_START_start << "\n";

#define DO_MEASURE_TIME_START time_t DO_MEASURE_TIME_START_start; time(&DO_MEASURE_TIME_START_start);
#define DO_MEASURE_TIME_END time_t DO_MEASURE_TIME_END_end; time(&DO_MEASURE_TIME_END_end); cout << "DO_MEASURE_TIME: " << DO_MEASURE_TIME_END_end - DO_MEASURE_TIME_START_start << "\n";


using namespace std;
using namespace cv;

typedef unsigned long ulong;
typedef unsigned int uint;

// return pair of index of min and max value
template<typename T>
pair<int,int> get_min_max(vector<T> vals) {
	if(vals.empty())
		throw new exception;

	int min_index = 0;
	T min = vals[0];

	int max_index = 0;
	T max = vals[0];

	for(int i=0; i<vals.size();i++) {
		const T &v = vals[i];
		if(v < min) {
			min = v;
			min_index = i;
		}
		if(v > max)
			max = v;
			max_index = i;
	}

	return make_pair(min_index,max_index);
}

void convert_glm_mat4_to_cv_Mat34(const glm::mat4 &glm_mat, cv::Mat &out);

void convert_glm_mat4_to_cv_Mat(const glm::mat4 &glm_mat, cv::Mat &out);

cv::Point2f convert_glm_vec2_to_cv_Point2f(const glm::vec2 &glm_vec);
cv::Point3f convert_glm_vec3_to_cv_Point3f(const glm::vec3 &glm_vec);

void convert_glm_vec2_to_cv_Point2f(const glm::vec2 &glm_vec, cv::Point2f &out);

std::string remove_path(std::string s);
std::string remove_extension(std::string s);

std::string remove_path_and_exception(std::string s);

void convert_glm_vec2_to_cv_Mat21(const glm::vec2 &glm_vec, cv::Mat &out);

vector<Mat> copy_imgs(const vector<Mat> &imgs);

vector<Mat> parse_images(const char *images_folder, const SfMData &sfm_data);
int parse_images(const char *images_folder, const SfMData &sfm_data, vector<Mat> &imgs);

vector<Mat> parse_imagesRGB(const char *images_folder, const SfMData &sfm_data);

vector<glm::vec4> convert_Vec4f_vec4(const vector<Vec4f> &cur_edges);

cv::Vec3f convert_from_glm_vec3_to_cv_Vec3f(const glm::vec3 &v, cv::Vec3f &vcv);

cv::Vec3f convert_from_glm_vec3_to_cv_Mat(const glm::vec3 &v, cv::Mat &vcv);

void convert_from_glm_mat3_to_cv_Mat3f(const glm::mat3 &m, cv::Mat &mcv);

set<int> find_points_on_both_images(const vector<set<int>> &points_on_images, int img_i, int img_j);

vector<pair<glm::vec2,glm::vec2>> find_fakepoints_on_both_images(const SfMData &sfmd, int img_i, int img_j, const ulong amount);

vector<set<int>> get_point_sets_on_images(const SfMData &sfmd);

vector<glm::vec2> get_vector_of_2d_coordinates_of_point_on_image(const SfMData &sfmd, const int img_id, const vector<int> &point_ids);

glm::vec2 get_2d_coordinates_of_point_on_image(const SfMData &sfmd, const int img_id, const int point_id, bool &found);

// Get color of 2D point on image using bilinear interpolation
cv::Vec3b getColorSubpix(const cv::Mat& img, cv::Point2f pt);

// Returns color of given point, by averaging color of all its 2D observations
Scalar get_refpoint_mixed_color(const SfMData &sfmd, const vector<Mat> &imgs, const ulong refpoint_id);

/**
 * Gets from SfMData the 2D coordinates of point_id on img_id
 */
void get_2d_coordinates_of_point_on_image(const SfMData &sfmd, const int img_id, const int point_id, glm::vec2 &point, bool &found);

template<typename T>
void print_Mat(cv::Mat m) {
	for(int i=0;i<m.rows;i++) {
		for(int j=0;j<m.cols;j++)
			cout << m.at<T>(i,j) << "\t";
		cout << endl;
	}
}
vector<Vec4f> convert_vec4_Vec4f(const vector<glm::vec4> &cur_edges);

int find_closest_point_on_image(const SfMData &sfmd, const int img_id, const glm::vec2 &p);

void print_vec2(const glm::vec2 &p2d);

void print_vec4(const glm::vec4 &p4d);

void print_vector_vec2(const vector<glm::vec2> &p2ds);

template<typename T>
void print_vector(vector<T> v) {
	cout << "[";
	if(v.size() > 0) {
		cout << v[0];
		for(int i=1;i<v.size();i++) {
			cout << ",";
			cout << v[i];
		}
	}
	cout << "]";
}

float avg_vec_float(const vector<float> &v);

void release_vector_Mat(vector<Mat> &matvec);

vector<Mat> copy_vector_Mat(vector<Mat> &matvec);

float avg_vec_float_weight_starting_image(const vector<float> &v, const vector<int> &selected_id, const int starting_img, const float starting_weight);

bool glm_vec2_equal(const glm::vec2 &v1, const glm::vec2 &v2);
bool vec_glm_vec2_equal(const vector<glm::vec2> &v1, const vector<glm::vec2> &v2);

bool glm_vec3_equal(const glm::vec3 &v1, const glm::vec3 &v2);
bool vec_glm_vec3_equal(const vector<glm::vec3> &v1, const vector<glm::vec3> &v2);

// check if two vector of vec2 are equal but reversed
bool vec_glm_vec2_equal_inv(const vector<glm::vec2> &v1, const vector<glm::vec2> &v2);
bool vec_glm_vec3_equal_inv(const vector<glm::vec3> &v1, const vector<glm::vec3> &v2);

std::ostream &operator<< (std::ostream &out, const glm::vec2 &vec);
std::ostream &operator<< (std::ostream &out, const glm::vec3 &vec);
std::ostream &operator<< (std::ostream &out, const glm::vec4 &vec);

float rand_float();

template<typename T>
std::ostream &operator<< (std::ostream &out, const vector<T> &vec)
{
	out << "[";
	if(vec.size() > 0) {
		out << vec[0];
		for(ulong i=1; i < vec.size(); i++)
			out << "," << vec[i];
	}
	out << "]";
	return out;
}

template<typename T1,typename T2>
std::ostream &operator<< (std::ostream &out, const pair<T1,T2> &s)
{
	return out << "(" << s.first << "," << s.second << ")";
}


template<typename T>
std::ostream &operator<< (std::ostream &out, const set<T> &s)
{
	out << "{";
	if(s.size() > 0) {
		int i=0;
		for(const auto &e : s)
			if(i==0) {
				out << e;
				i++;
			} else
				out << "," << e;
	}
	out << "}";
	return out;
}

std::vector<glm::vec2> intersections_remove_segments(const std::vector<SegmentEdgeManager::segment_point> &intersections);

std::vector<std::vector<std::vector<glm::vec2>>> correspondences_remove_segments(const std::vector<std::vector<std::vector<SegmentEdgeManager::segment_point>>> &correspondences);

template<typename T>
bool is_in(const set<T> &s, const T &e) {
	return s.find(e) != s.end();
}

template<typename A, typename B>
vector<pair<A,B>> zip_pair(const vector<A> &a, const vector<B> &b) {
	vector<pair<A,B>> res;
	for(int i=0; i < min(a.size(),b.size()); i++)
		res.push_back(make_pair(a[i],b[i]));
	return res;
}

// reorder A and B according to B's values
template<typename A, typename B>
void reorder_pair_of_vector(vector<A> &a, vector<B> &b) {
	vector<pair<A,B>> vab = zip_pair(a,b);

	struct ordering {
	    bool operator ()(pair<A, B> const& a,
	                     pair<A, B> const& b) {
	        return a.second < b.second;
	    }
	};

	std::sort(vab.begin(),vab.end(),ordering());

	for(int i=0; i < vab.size(); i++) {
		a[i] = vab[i].first;
		b[i] = vab[i].second;
	}
}

glm::mat4x4 get_rt4x4(const glm::mat3 &r, const glm::vec3 &t);

template<typename A, typename B, typename C>
std::ostream &operator<< (std::ostream &out, const tuple<A,B,C> abc) {
	out << "(" << get<0>(abc) << "," << get<1>(abc) << "," << get<2>(abc) << ")";
	return out;
}

template<typename T>
T vec_max(vector<T> v) {
	if(v.size() == 0)
		std::invalid_argument("vec_max : zero-sized vector!");
	T maxe = v[0];
	for(const auto &e : v)
		maxe = maxe >= e ? maxe : e;

	return maxe;
}

std::ostream &operator<<(std::ostream &out, const glm::mat3 &m);
std::ostream &operator<<(std::ostream &out, const glm::mat4 &m);

template<typename T>
T** create_2D_array(unsigned long nrows, unsigned long ncols) {
	T** ptr = new T*[nrows];  // allocate pointers
	T* pool = new T[nrows * ncols];  // allocate pool
	for (unsigned i = 0; i < nrows; ++i, pool += ncols)
		ptr[i] = pool;
	return ptr;
}

template<typename T>
void delete_2D_array(T** arr) {
	delete[] arr[0];  // remove the pool
	delete[] arr;     // remove the pointers
}

Mat** create_2D_Mat_array(unsigned long nrows, unsigned long ncols);

void delete_2D_Mat_array(Mat** arr, unsigned long nrows, unsigned long ncols);

float floor_or_upper_if_close(const float v);
bool is_m_multiple_of_n_float(const float m, const float n);

// valid is set to false if the coords are on a boundary
pair<ulong,ulong> get_2dmap_cell_from_coords(const float cell_dim, const glm::vec2 &coords, bool &on_boundary);
pair<ulong,ulong> get_2dmap_cell_from_coords(const float cell_dim, const glm::vec2 &coords, bool &on_boundary_row, bool &on_boundary_col);

template<typename T1>
cv::Vec<T1, 3> convert_vec3_Vec3(glm::vec3 &v) {
	return cv::Vec<T1, 3>(v[0],v[1],v[2]);
}

template<typename T1>
vector<cv::Vec<T1, 3>> convert_vectorvec3_vectorVec3(const vector<glm::vec3> &v) {
	vector<cv::Vec<T1, 3>> vV;
	for(const auto &ve: v)
		vV.push_back(convert_vec3_Vec3<T1>(ve));
	return vV;
}

string compute_plg_path(char *out_images_folder, const SfMData &sfm_data, int img_id, const String &s);

string compute_img_path(char *out_images_folder, const SfMData &sfm_data, int img_id, const String &s);

void write_img(char *out_images_folder, const SfMData &sfm_data, int img_id, const Mat &img, const String &s);

int write_images(const string &out_images_folder, const SfMData &sfm_data, vector<Mat> &imgs, const String &s);

int write_images(const string &out_images_folder, const SfMData &sfm_data, vector<Mat> &imgs);

pair<vector<glm::vec3>,vector<Scalar>> parseGTxyz(const string &spath);

pair<vector<glm::vec3>,vector<Scalar>> parseGTply(const string &spath);

glm::vec3 to_glmvec3(CGAL::Simple_cartesian<double>::Point_3 p);

std::ostream &operator<< (std::ostream &out, const pair<vector<int>,vector<ulong>> &p);

SfMData read_sfm_data(const char *sfmd_path);

#endif /* INCLUDE_EDGEGRAPH3D_UTILS_EDGE_GRAPH_3D_UTILITIES_HPP_ */
