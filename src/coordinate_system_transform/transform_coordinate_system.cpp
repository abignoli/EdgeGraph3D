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

#include "transform_coordinate_system.hpp"

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>

#include "SfMData.h"
#include "types_reconstructor.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "glm.hpp"


using namespace Eigen;
using namespace std;

// given sfmd and file containing actual camera positions, update sfmd
void test_umeyama() {
	Vector3d a1(1,1,0);
	Vector3d a2(-1,1,0);
	Vector3d a3(-1,-1,0);
	Vector3d a4(1,-1,0);

	Vector3d b1(2,2,0);
	Vector3d b2(-2,2,0);
	Vector3d b3(-2,-2,0);
	Vector3d b4(2,-2,0);

	Matrix<double,3,4> start,end;
	start.col(0)=a1;
	start.col(1)=a2;
	start.col(2)=a3;
	start.col(3)=a4;


	end.col(0)=b1;
	end.col(1)=b2;
	end.col(2)=b3;
	end.col(3)=b4;

	cout << Eigen::umeyama(start,end,true) << endl;
}

MatrixXd read_camera_file(const std::string &cameras_file, const int numCams) {
	MatrixXd cam_pos(3,numCams);
	std::ifstream myf;
	myf.open(cameras_file);
	for(int i=0; i < numCams;i++) {
		//std::string sparam;
		double a,b,c;
		myf >> a;
		myf >> b;
		myf >> c;

		cam_pos.col(i)=Vector3d(a,b,c);
	}
	myf.close();
	return cam_pos;
}

Matrix<double,4,4> compute_transformation(SfMData &sfmd, const std::string &cameras_file) {
	MatrixXd real_cam_pos = read_camera_file(cameras_file,sfmd.numCameras_);
	MatrixXd cur_cam_pos(3,sfmd.numCameras_);
	for(int i=0; i < sfmd.numCameras_; i++)
		cur_cam_pos.col(i)=Vector3d(sfmd.camerasList_[i].center.x,sfmd.camerasList_[i].center.y,sfmd.camerasList_[i].center.z);
	return Eigen::umeyama(cur_cam_pos,real_cam_pos,true);
}

bool is_nullCamera(const CameraType &ct) {
	return  (ct.center.x==0.0) && (ct.center.y==0.0) && (ct.center.z==0.0) &&
			(ct.rotation[0][0]==0.0) && (ct.rotation[0][1]==0.0) && (ct.rotation[0][2]==0.0) &&
			(ct.rotation[1][0]==0.0) && (ct.rotation[1][1]==0.0) && (ct.rotation[1][2]==0.0) &&
			(ct.rotation[2][0]==0.0) && (ct.rotation[2][1]==0.0) && (ct.rotation[2][2]==0.0);
}

set<int> compute_nullCameras(SfMData &sfmd) {
	set<int> res;
	for(int i=0; i < sfmd.numCameras_; i++)
		if(is_nullCamera(sfmd.camerasList_[i]))
			res.insert(i);
	return res;
}

MatrixXd read_camera_file_exclude_nullCamera(const std::string &cameras_file, const int numCams, const set<int> &nullCameras) {
	MatrixXd cam_pos(3,numCams- nullCameras.size());
	std::ifstream myf;
	myf.open(cameras_file);
	int skipped=0;
	for(int i=0; i < numCams;i++) {
		double a,b,c;
		myf >> a;
		myf >> b;
		myf >> c;
		if(nullCameras.find(i) == nullCameras.end()) {
			//std::string sparam;
			cam_pos.col(i-skipped)=Vector3d(a,b,c);
		} else
			skipped++;
	}

	myf.close();
	return cam_pos;
}

Matrix<double,4,4> compute_transformation_exclude_nullCamera(SfMData &sfmd, const std::string &cameras_file) {
	set<int> nullCameras = compute_nullCameras(sfmd);
	MatrixXd real_cam_pos = read_camera_file_exclude_nullCamera(cameras_file,sfmd.numCameras_,nullCameras);
	MatrixXd cur_cam_pos(3,sfmd.numCameras_- nullCameras.size());
	int skipped=0;
	for(int i=0; i < sfmd.numCameras_; i++)
		if(nullCameras.find(i) == nullCameras.end())
			cur_cam_pos.col(i-skipped)=Vector3d(sfmd.camerasList_[i].center.x,sfmd.camerasList_[i].center.y,sfmd.camerasList_[i].center.z);
		else {
			cout << "Ignoring camera " << i << " in computing the Umeyama transformation" << endl;
			skipped++;
		}
	return Eigen::umeyama(cur_cam_pos,real_cam_pos,true);
}

glm::vec3 transform_3d_point(const glm::vec3 &pt, const glm::mat4 &transformation) {
	glm::vec4 tmp = glm::vec4(pt,1)*transformation;
	return glm::vec3(tmp[0]/tmp[3],tmp[1]/tmp[3],tmp[2]/tmp[3]);
}

void update_camera_transform(CameraType &c, const glm::mat4 &transformation) {
	if(is_nullCamera(c))
		return;

	c.cameraMatrix *= transformation;

	glm::mat4 e = glm::inverse(transformation) * get_rt4x4(c.rotation,c.translation);
	glm::mat4 kMatrix(0.0);
    for (int curR = 0; curR < 3; curR++) {
      for (int curC = 0; curC < 3; curC++) {
        kMatrix[curR][curC] = c.intrinsics[curR][curC];
      }
    }
	c.cameraMatrix = e * kMatrix;
	c.rotation = glm::mat3(e[0][0],e[0][1],e[0][2],e[1][0],e[1][1],e[1][2],e[2][0],e[2][1],e[2][2]);
	c.translation = glm::vec3(e[0][3],e[1][3],e[2][3]);
	c.center = transform_3d_point(c.center, transformation);
}

float compute_camera_errors(SfMData &sfmd, const std::string &cameras_file) {
	MatrixXd real_cam_pos = read_camera_file(cameras_file,sfmd.numCameras_);
	MatrixXd cur_cam_pos(3,sfmd.numCameras_);
	for(int i=0; i < sfmd.numCameras_; i++)
		cur_cam_pos.col(i)=Vector3d(sfmd.camerasList_[i].center.x,sfmd.camerasList_[i].center.y,sfmd.camerasList_[i].center.z);

	MatrixXd diffs = cur_cam_pos-real_cam_pos;
	float sum=0;
	for(int i=0; i < sfmd.numCameras_; i++)
		sum+=diffs.col(i).norm();
	return sum / sfmd.numCameras_;
}

float compute_camera_errors_exclude_nullCamera(SfMData &sfmd, const std::string &cameras_file) {
	set<int> nullCameras = compute_nullCameras(sfmd);
	MatrixXd real_cam_pos = read_camera_file_exclude_nullCamera(cameras_file,sfmd.numCameras_,nullCameras);
	MatrixXd cur_cam_pos(3,sfmd.numCameras_- nullCameras.size());
	int skipped=0;
	for(int i=0; i < sfmd.numCameras_; i++)
		if(nullCameras.find(i) == nullCameras.end())
			cur_cam_pos.col(i-skipped)=Vector3d(sfmd.camerasList_[i].center.x,sfmd.camerasList_[i].center.y,sfmd.camerasList_[i].center.z);
		else {
			cout << "Ignoring camera " << i << "in computing the Umeyama transformation" << endl;
			skipped++;
		}
	MatrixXd diffs = cur_cam_pos-real_cam_pos;
	float sum=0;
	for(int i=0; i < sfmd.numCameras_- nullCameras.size(); i++)
		sum+=diffs.col(i).norm();
	return sum / sfmd.numCameras_;
}

void update_sfmd_transform(SfMData &sfmd, const glm::mat4 &transformation) {
	// transform cameras
	for(int i=0;i<sfmd.numCameras_;i++)
		update_camera_transform(sfmd.camerasList_[i],transformation);

	// transform points
	for(int refpoint=0;refpoint<sfmd.numPoints_;refpoint++)
		sfmd.points_[refpoint]=transform_3d_point(sfmd.points_[refpoint],transformation);
}



void update_sfmd_from_real_camera_positions_no_null_cam_handle(SfMData &sfmd, const std::string &cameras_file) {
	Matrix<double,4,4> transformation = compute_transformation(sfmd,cameras_file);

	glm::mat4 glm_transformation;
	for(int i=0; i < 4; i++)
		for(int j=0; j < 4; j++)
			glm_transformation[i][j] = transformation(i,j);

	cout << "Camera position errors before: " << compute_camera_errors(sfmd,cameras_file) << endl;
	update_sfmd_transform(sfmd,glm_transformation);
	cout << "Camera position errors after: " << compute_camera_errors(sfmd,cameras_file) << endl;

	//cout << transformation;
}



void update_sfmd_from_real_camera_positions(SfMData &sfmd, const std::string &cameras_file) {

	cout << "Transforming coordinate system to adapt it to real camera poses..." << endl;
	Matrix<double,4,4> transformation = compute_transformation_exclude_nullCamera(sfmd,cameras_file);

	cout << "T:" << endl;
	cout << transformation << endl;;

	glm::mat4 glm_transformation;
	for(int i=0; i < 4; i++)
		for(int j=0; j < 4; j++)
			glm_transformation[i][j] = transformation(i,j);

	cout << "Camera position errors before: " << compute_camera_errors_exclude_nullCamera(sfmd,cameras_file) << endl;
	update_sfmd_transform(sfmd,glm_transformation);
	cout << "Camera position errors after: " << compute_camera_errors_exclude_nullCamera(sfmd,cameras_file) << endl;

	//cout << transformation;
}
