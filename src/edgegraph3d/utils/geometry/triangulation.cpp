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


#include "triangulation.hpp"

#include <opencv2/calib3d.hpp>
#include <opencv2/core/cvdef.h>
#include <opencv2/core/mat.hpp>
#include <opencv2/core/mat.inl.hpp>
#include <opencv2/core/matx.hpp>
#include <opencv2/core/types.hpp>
#include <opencv2/core.hpp>
#include <algorithm>
#include <cmath>
#include <utility>

#include <glm.hpp>
#include "SfMData.h"
#include "types_reconstructor.hpp"
#include "plg_matching.hpp"
#include "edge_graph_3d_utilities.hpp"


using namespace std;

int em_point2D3DJacobian(const std::vector<cv::Mat> &cameras, const cv::Mat &cur3Dpoint, cv::Mat &J, cv::Mat &hessian) {

  int numMeasures = cameras.size();
  cv::Mat cur3DPointHomog = cv::Mat(4, 1, CV_64F);

  cur3DPointHomog.at<double>(0, 0) = cur3Dpoint.at<double>(0, 0);
  cur3DPointHomog.at<double>(1, 0) = cur3Dpoint.at<double>(1, 0);
  cur3DPointHomog.at<double>(2, 0) = cur3Dpoint.at<double>(2, 0);
  cur3DPointHomog.at<double>(3, 0) = 1.0;

  J = cv::Mat(2 * numMeasures, 3, CV_64FC1);  //2 rows for each point: one for x, the other for y
  hessian = cv::Mat(3, 3, CV_64FC1);

  for (int curMeas = 0; curMeas < numMeasures; ++curMeas) {
    cv::Mat curReproj = cameras[curMeas] * cur3DPointHomog;
    double xH = curReproj.at<double>(0, 0);
    double yH = curReproj.at<double>(1, 0);
    double zH = curReproj.at<double>(2, 0);
    double p00 = cameras[curMeas].at<double>(0, 0);
    double p01 = cameras[curMeas].at<double>(0, 1);
    double p02 = cameras[curMeas].at<double>(0, 2);
    double p10 = cameras[curMeas].at<double>(1, 0);
    double p11 = cameras[curMeas].at<double>(1, 1);
    double p12 = cameras[curMeas].at<double>(1, 2);
    double p20 = cameras[curMeas].at<double>(2, 0);
    double p21 = cameras[curMeas].at<double>(2, 1);
    double p22 = cameras[curMeas].at<double>(2, 2);

    //d(P*X3D)/dX
    J.at<double>(2 * curMeas, 0) = (p00 * zH - p20 * xH) / (zH * zH);
    J.at<double>(2 * curMeas + 1, 0) = (p10 * zH - p20 * yH) / (zH * zH);

    //d(P*X3D)/dY
    J.at<double>(2 * curMeas, 1) = (p01 * zH - p21 * xH) / (zH * zH);
    J.at<double>(2 * curMeas + 1, 1) = (p11 * zH - p21 * yH) / (zH * zH);

    //d(P*X3D)/dZ
    J.at<double>(2 * curMeas, 2) = (p02 * zH - p22 * xH) / (zH * zH);
    J.at<double>(2 * curMeas + 1, 2) = (p12 * zH - p22 * yH) / (zH * zH);
  }

  hessian = J.t() * J;
  double d;
  d = cv::determinant(hessian);
  if (d < 0.00001) {
    //// printf("doh");
    return -1;
  } else {
    return 1;
  }
}

int em_GaussNewton(const std::vector<cv::Mat> &cameras, const std::vector<cv::Point2f> &points, const cv::Point3d &init3Dpoint,
  cv::Point3d &optimizedPoint) {

  int numMeasures = points.size();
  cv::Mat r = cv::Mat(numMeasures * 2, 1, CV_64F);
  cv::Mat curEstimate3DPoint = cv::Mat(3, 1, CV_64F);
  cv::Mat curEstimate3DPointH = cv::Mat(4, 1, CV_64F);

  curEstimate3DPoint.at<double>(0, 0) = init3Dpoint.x;
  curEstimate3DPoint.at<double>(1, 0) = init3Dpoint.y;
  curEstimate3DPoint.at<double>(2, 0) = init3Dpoint.z;

  cv::Mat J, H;
  double last_mse = 0;

  /************NB maxGaussNewtonIteration needs to be fixed ***********/
  // number of iterations set to 30
  for (int i = 0; i < 30; i++) {

    double mse = 0;
        //compute residuals
    for (int curMeas = 0; curMeas < numMeasures; ++curMeas) {
      curEstimate3DPointH.at<double>(0, 0) = curEstimate3DPoint.at<double>(0, 0);
      curEstimate3DPointH.at<double>(1, 0) = curEstimate3DPoint.at<double>(1, 0);
      curEstimate3DPointH.at<double>(2, 0) = curEstimate3DPoint.at<double>(2, 0);
      curEstimate3DPointH.at<double>(3, 0) = 1.0;
      cv::Mat cur2DpositionH = cameras[curMeas] * curEstimate3DPointH;

      r.at<double>(2 * curMeas, 0) = ((points[curMeas].x - cur2DpositionH.at<double>(0, 0) / cur2DpositionH.at<double>(2, 0)));
      mse += r.at<double>(2 * curMeas, 0) * r.at<double>(2 * curMeas, 0);

      r.at<double>(2 * curMeas + 1, 0) = ((points[curMeas].y - cur2DpositionH.at<double>(1, 0) / cur2DpositionH.at<double>(2, 0)));
      mse += r.at<double>(2 * curMeas + 1, 0) * r.at<double>(2 * curMeas + 1, 0);
#ifdef DEBUG_OPTIMIZATION_VERBOSE
      std::cout<<"CurMeas: "<<curMeas<<std::endl<<"curEstimate3DPointH="<< curEstimate3DPointH.t()<<std::endl;
      std::cout<<"CurCam"<<cameras[curMeas]<<std::endl;
      std::cout<<"cur2DpositionH: "<<cur2DpositionH.at<double>(0, 0)/cur2DpositionH.at<double>(2, 0)<<", "<<cur2DpositionH.at<double>(1, 0)/cur2DpositionH.at<double>(2, 0)<<std::endl;
      std::cout<<"points[curMeas]: "<<points[curMeas]<<std::endl;
      std::cout<<"residual on x: "<<r.at<double>(2 * curMeas, 0)<<std::endl;
      std::cout<<"residual on y: "<<r.at<double>(2 * curMeas + 1 , 0)<<std::endl;
      std::cout<<std::endl;
#endif
    }

// if the error is very low, it  ends the function
    if (abs(mse / (numMeasures * 2) - last_mse) < 0.0000005) {
      break;
    }
    last_mse = mse / (numMeasures * 2);

    if (em_point2D3DJacobian(cameras, curEstimate3DPoint, J, H) == -1)
      return -1;
#ifdef DEBUG_OPTIMIZATION_VERBOSE
    std::cout<<"J: "<<J<<std::endl;
    std::cout<<"H: "<<H<<std::endl;
#endif

    curEstimate3DPoint += H.inv() * J.t() * r;

#ifdef DEBUG_OPTIMIZATION
    // printf("%d %f\n", i, last_mse);
#endif
  }
  if (last_mse < 9/*3 pixels*/) {
    optimizedPoint.x = curEstimate3DPoint.at<double>(0, 0);
    optimizedPoint.y = curEstimate3DPoint.at<double>(1, 0);
    optimizedPoint.z = curEstimate3DPoint.at<double>(2, 0);
    return 1;
  } else {
    return -1;
  }
}

void em_estimate3Dpositions(const SfMData &sfmd, const vector<glm::vec2> &selected_2d_reprojections_coords, const vector<int> &selected_2d_reprojections_ids, glm::vec3 &triangulated_point, bool &valid) {
	const pair<int, int> min_max_ind = get_min_max(
			selected_2d_reprojections_ids);

	/*
	 * Note: I decided to use the minimum and maximum indices for the selection
	 * of the first and last cameras used for initial triangulation. This is done with the
	 * assumption that the camera rays between those cameras form a bigger
	 * angle than the other possible camera pairs. The reason for this is that
	 * I suspect that having wider angle between the cameras used for initialization,
	 * will lead to a better initial result.
	 */

	int firstCamIdx, lastCamIdx;
	cv::Point3d triangulated3DPointInit, triangulated3DPoint;
	cv::Vec4f triangulated3DPointInitTemp;
	vector<cv::Point2f> firstPositionVec, lastPositionVec;
	cv::Mat firstCam, lastCam;
	cv::Point2f firstPoint, lastPoint;
	std::vector<cv::Mat> curCams;
	std::vector<cv::Point2f> curPoints;

	firstCamIdx = selected_2d_reprojections_ids[min_max_ind.first];
	lastCamIdx = selected_2d_reprojections_ids[min_max_ind.second];

	convert_glm_mat4_to_cv_Mat34(sfmd.camerasList_[firstCamIdx].cameraMatrix,
			firstCam);
	convert_glm_mat4_to_cv_Mat34(sfmd.camerasList_[lastCamIdx].cameraMatrix,
			lastCam);

	convert_glm_vec2_to_cv_Point2f(
			selected_2d_reprojections_coords[min_max_ind.first], firstPoint);
	convert_glm_vec2_to_cv_Point2f(selected_2d_reprojections_coords[min_max_ind.second],
			lastPoint);

	firstPositionVec.push_back(firstPoint);
	lastPositionVec.push_back(lastPoint);

	cv::triangulatePoints(firstCam, lastCam, firstPositionVec, lastPositionVec,
			triangulated3DPointInitTemp);

	triangulated3DPointInit.x = triangulated3DPointInitTemp[0]
			/ triangulated3DPointInitTemp[3];
	triangulated3DPointInit.y = triangulated3DPointInitTemp[1]
			/ triangulated3DPointInitTemp[3];
	triangulated3DPointInit.z = triangulated3DPointInitTemp[2]
			/ triangulated3DPointInitTemp[3];

	//Pack the information to start the Gauss Newton algorithm
	for (int i = 0; i < selected_2d_reprojections_ids.size(); i++) {
		// FIXME change this creating directly CV_64F
		int camIdx = selected_2d_reprojections_ids[i];
		cv::Mat temp1, temp2;
		convert_glm_mat4_to_cv_Mat(sfmd.camerasList_[camIdx].cameraMatrix, temp1);
		temp1.convertTo(temp2,CV_64F);
		curCams.push_back(temp2);

		cv::Point2f curPoint;
		convert_glm_vec2_to_cv_Point2f(
					selected_2d_reprojections_coords[i], curPoint);
		curPoints.push_back(curPoint);
	}

	int resGN = em_GaussNewton(curCams, curPoints, triangulated3DPointInit, triangulated3DPoint);

	if(resGN != -1) {
		triangulated_point[0] = triangulated3DPoint.x;
		triangulated_point[1] = triangulated3DPoint.y;
		triangulated_point[2] = triangulated3DPoint.z;
		valid = true;
	} else
		valid = false;
}

void em_estimate3Dpositions(const SfMData &sfmd, const vector<PolyLineGraph2D::plg_point> &selected_2d_reprojections_coords, const vector<int> &selected_2d_reprojections_ids, glm::vec3 &triangulated_point, bool &valid) {
	const pair<int, int> min_max_ind = get_min_max(
			selected_2d_reprojections_ids);

	/*
	 * Note: I decided to use the minimum and maximum indices for the selection
	 * of the first and last cameras used for initial triangulation. This is done with the
	 * assumption that the camera rays between those cameras form a bigger
	 * angle than the other possible camera pairs. The reason for this is that
	 * I suspect that having wider angle between the cameras used for initialization,
	 * will lead to a better initial result.
	 */

	int firstCamIdx, lastCamIdx;
	cv::Point3d triangulated3DPointInit, triangulated3DPoint;
	cv::Vec4f triangulated3DPointInitTemp;
	vector<cv::Point2f> firstPositionVec, lastPositionVec;
	cv::Mat firstCam, lastCam;
	cv::Point2f firstPoint, lastPoint;
	std::vector<cv::Mat> curCams;
	std::vector<cv::Point2f> curPoints;

	firstCamIdx = selected_2d_reprojections_ids[min_max_ind.first];
	lastCamIdx = selected_2d_reprojections_ids[min_max_ind.second];

	convert_glm_mat4_to_cv_Mat34(sfmd.camerasList_[firstCamIdx].cameraMatrix,
			firstCam);
	convert_glm_mat4_to_cv_Mat34(sfmd.camerasList_[lastCamIdx].cameraMatrix,
			lastCam);

	convert_glm_vec2_to_cv_Point2f(
			selected_2d_reprojections_coords[min_max_ind.first].plp.coords, firstPoint);
	convert_glm_vec2_to_cv_Point2f(selected_2d_reprojections_coords[min_max_ind.second].plp.coords,
			lastPoint);

	firstPositionVec.push_back(firstPoint);
	lastPositionVec.push_back(lastPoint);

	cv::triangulatePoints(firstCam, lastCam, firstPositionVec, lastPositionVec,
			triangulated3DPointInitTemp);

	triangulated3DPointInit.x = triangulated3DPointInitTemp[0]
			/ triangulated3DPointInitTemp[3];
	triangulated3DPointInit.y = triangulated3DPointInitTemp[1]
			/ triangulated3DPointInitTemp[3];
	triangulated3DPointInit.z = triangulated3DPointInitTemp[2]
			/ triangulated3DPointInitTemp[3];

	//Pack the information to start the Gauss Newton algorithm
	for (int i = 0; i < selected_2d_reprojections_ids.size(); i++) {
		int camIdx = selected_2d_reprojections_ids[i];
		cv::Mat temp1, temp2;
		convert_glm_mat4_to_cv_Mat(sfmd.camerasList_[camIdx].cameraMatrix, temp1);
		temp1.convertTo(temp2,CV_64F);
		curCams.push_back(temp2);

		cv::Point2f curPoint;
		convert_glm_vec2_to_cv_Point2f(
					selected_2d_reprojections_coords[i].plp.coords, curPoint);
		curPoints.push_back(curPoint);
	}

	int resGN = em_GaussNewton(curCams, curPoints, triangulated3DPointInit, triangulated3DPoint);

	if(resGN != -1) {
		triangulated_point[0] = triangulated3DPoint.x;
		triangulated_point[1] = triangulated3DPoint.y;
		triangulated_point[2] = triangulated3DPoint.z;
		valid = true;
	} else
		valid = false;
}

bool compatible_new_observation_to_3Dpositions(const SfMData &sfmd, const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &current_point, const glm::vec2 new_coords, const int new_viewpoint_id, glm::vec3 &triangulated_point) {
	bool valid;
	em_add_new_observation_to_3Dpositions(sfmd, current_point, new_coords, new_viewpoint_id, triangulated_point, valid);
	return valid;
}

bool compatible_new_observation_to_3Dpositions(const SfMData &sfmd, const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &current_point, const PolyLineGraph2D::plg_point &new_plgp, const int new_viewpoint_id, glm::vec3 &triangulated_point) {
	return compatible_new_observation_to_3Dpositions(sfmd, current_point, new_plgp.plp.coords, new_viewpoint_id, triangulated_point);
}

bool compatible_new_observation_to_3Dpositions_update(const SfMData &sfmd, std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &current_point, const PolyLineGraph2D::plg_point &new_plgp, const int new_viewpoint_id) {
	glm::vec3 triangulated_point;
	bool valid = compatible_new_observation_to_3Dpositions(sfmd, current_point, new_plgp, new_viewpoint_id, triangulated_point);
	if(!valid)
		return false;
	// update point
	get<0>(current_point) = triangulated_point;
	get<1>(current_point).push_back(new_plgp);
	get<2>(current_point).push_back(new_viewpoint_id);
	return true;
}

void em_add_new_observation_to_3Dpositions(const SfMData &sfmd, const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &current_point, const glm::vec2 new_coords, const int new_viewpoint_id, glm::vec3 &triangulated_point, bool &valid) {


	/*
	 * Note: I decided to use the minimum and maximum indices for the selection
	 * of the first and last cameras used for initial triangulation. This is done with the
	 * assumption that the camera rays between those cameras form a bigger
	 * angle than the other possible camera pairs. The reason for this is that
	 * I suspect that having wider angle between the cameras used for initialization,
	 * will lead to a better initial result.
	 */

	cv::Point3d triangulated3DPointInit, triangulated3DPoint;
	cv::Vec4f triangulated3DPointInitTemp;
	vector<cv::Point2f> firstPositionVec, lastPositionVec;
	cv::Mat firstCam, lastCam;
	cv::Point2f firstPoint, lastPoint;
	std::vector<cv::Mat> curCams;
	std::vector<cv::Point2f> curPoints;

	triangulated3DPointInit.x = get<0>(current_point)[0];
	triangulated3DPointInit.y = get<0>(current_point)[1];
	triangulated3DPointInit.z = get<0>(current_point)[2];

	//Pack the information to start the Gauss Newton algorithm
	for (int i = 0; i < get<2>(current_point).size(); i++) {
		// FIXME change this creating directly CV_64F
		int camIdx = get<2>(current_point)[i];
		cv::Mat temp1, temp2;
		convert_glm_mat4_to_cv_Mat(sfmd.camerasList_[camIdx].cameraMatrix, temp1);
		temp1.convertTo(temp2,CV_64F);
		curCams.push_back(temp2);

		cv::Point2f curPoint;
		convert_glm_vec2_to_cv_Point2f(
				get<1>(current_point)[i].plp.coords, curPoint);
		curPoints.push_back(curPoint);
	}

	int camIdx = new_viewpoint_id;
	cv::Mat temp1, temp2;
	convert_glm_mat4_to_cv_Mat(sfmd.camerasList_[camIdx].cameraMatrix, temp1);
	temp1.convertTo(temp2,CV_64F);
	curCams.push_back(temp2);

	cv::Point2f curPoint;
	convert_glm_vec2_to_cv_Point2f(new_coords,curPoint);
	curPoints.push_back(curPoint);

	int resGN = em_GaussNewton(curCams, curPoints, triangulated3DPointInit, triangulated3DPoint);

	if(resGN != -1) {
		triangulated_point[0] = triangulated3DPoint.x;
		triangulated_point[1] = triangulated3DPoint.y;
		triangulated_point[2] = triangulated3DPoint.z;
		valid = true;
	} else
		valid = false;
}


void em_add_new_observation_to_3Dpositions(const SfMData &sfmd, const glm::vec3 &current_point_coords, const vector<glm::vec2> &current_point_observation_coords, const vector<int> &current_point_observation_ids, const glm::vec2 new_coords, const int new_viewpoint_id, glm::vec3 &triangulated_point, bool &valid) {


	/*
	 * Note: I decided to use the minimum and maximum indices for the selection
	 * of the first and last cameras used for initial triangulation. This is done with the
	 * assumption that the camera rays between those cameras form a bigger
	 * angle than the other possible camera pairs. The reason for this is that
	 * I suspect that having wider angle between the cameras used for initialization,
	 * will lead to a better initial result.
	 */

	cv::Point3d triangulated3DPointInit, triangulated3DPoint;
	cv::Vec4f triangulated3DPointInitTemp;
	vector<cv::Point2f> firstPositionVec, lastPositionVec;
	cv::Mat firstCam, lastCam;
	cv::Point2f firstPoint, lastPoint;
	std::vector<cv::Mat> curCams;
	std::vector<cv::Point2f> curPoints;

	triangulated3DPointInit.x = current_point_coords[0];
	triangulated3DPointInit.y = current_point_coords[1];
	triangulated3DPointInit.z = current_point_coords[2];

	//Pack the information to start the Gauss Newton algorithm
	for (int i = 0; i < current_point_observation_ids.size(); i++) {
		// FIXME change this creating directly CV_64F
		int camIdx = current_point_observation_ids[i];
		cv::Mat temp1, temp2;
		convert_glm_mat4_to_cv_Mat(sfmd.camerasList_[camIdx].cameraMatrix, temp1);
		temp1.convertTo(temp2,CV_64F);
		curCams.push_back(temp2);

		cv::Point2f curPoint;
		convert_glm_vec2_to_cv_Point2f(
				current_point_observation_coords[i], curPoint);
		curPoints.push_back(curPoint);
	}

	int camIdx = new_viewpoint_id;
	cv::Mat temp1, temp2;
	convert_glm_mat4_to_cv_Mat(sfmd.camerasList_[camIdx].cameraMatrix, temp1);
	temp1.convertTo(temp2,CV_64F);
	curCams.push_back(temp2);

	cv::Point2f curPoint;
	convert_glm_vec2_to_cv_Point2f(new_coords,curPoint);
	curPoints.push_back(curPoint);

	int resGN = em_GaussNewton(curCams, curPoints, triangulated3DPointInit, triangulated3DPoint);

	if(resGN != -1) {
		triangulated_point[0] = triangulated3DPoint.x;
		triangulated_point[1] = triangulated3DPoint.y;
		triangulated_point[2] = triangulated3DPoint.z;
		valid = true;
	} else
		valid = false;
}

void compute_3d_point_coords(const SfMData &sfmd,
		const vector<glm::vec2> &selected_2d_reprojections_coords,
		const vector<int> &selected_2d_reprojections_ids,
		glm::vec3 &new_point_data,
		bool &valid) {
	valid = false;

	if (selected_2d_reprojections_coords.size() >= 2) {
		// Compute 3D point
		em_estimate3Dpositions(sfmd, selected_2d_reprojections_coords,
				selected_2d_reprojections_ids, new_point_data,
				valid);
	}
}

void compute_3d_point(const SfMData &sfmd,
		const vector<PolyLineGraph2D::plg_point> &selected_2d_reprojections_coords,
		const vector<int> &selected_2d_reprojections_ids,
		std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &new_point_data,
		bool &valid) {
	valid = false;

	if (selected_2d_reprojections_coords.size() >= 2) {
		em_estimate3Dpositions(sfmd, selected_2d_reprojections_coords,
				selected_2d_reprojections_ids, std::get < 0 > (new_point_data),
				valid);
		if (valid) {
			std::get < 1 > (new_point_data) = selected_2d_reprojections_coords;
			std::get < 2 > (new_point_data) = selected_2d_reprojections_ids;
		}
	}
}

void compute_all_potential_3d_points_3views(const SfMData &sfmd,
		const tuple<vector<PolyLineGraph2D::plg_point>,vector<PolyLineGraph2D::plg_point>,vector<PolyLineGraph2D::plg_point>> &all_2d_reprojections_3views,
		const tuple<int,int,int> &all_2d_reprojections_ids_3views,
		vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &potential_new_points_data) {
	potential_new_points_data.clear();

	bool valid;
	vector<int> views;
	views.push_back(get<0>(all_2d_reprojections_ids_3views));
	views.push_back(get<1>(all_2d_reprojections_ids_3views));
	views.push_back(get<2>(all_2d_reprojections_ids_3views));

	for(const auto &plgp0 : get<0>(all_2d_reprojections_3views))
		for(const auto &plgp1 : get<1>(all_2d_reprojections_3views))
			for(const auto &plgp2 : get<2>(all_2d_reprojections_3views)) {
				glm::vec3 triangulated_point;
				vector<PolyLineGraph2D::plg_point> new_plgps;
				new_plgps.push_back(plgp0);
				new_plgps.push_back(plgp1);
				new_plgps.push_back(plgp2);

				em_estimate3Dpositions(sfmd, new_plgps, views, triangulated_point, valid);

				vector<int> views_cpy = views;
				if(valid)
					potential_new_points_data.push_back(make_tuple(triangulated_point,new_plgps,views_cpy));
			}
}

void compute_all_potential_3d_points_3views_plg_following_newpoint_compatibility(const SfMData &sfmd,const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs,
		const tuple<vector<PolyLineGraph2D::plg_point>,vector<PolyLineGraph2D::plg_point>,vector<PolyLineGraph2D::plg_point>> &all_2d_reprojections_3views,
		const tuple<int,int,int> &all_2d_reprojections_ids_3views,
		vector<new_3dpoint_and_sides_plgp_matches> &potential_new_points_and_following_data_directions) {
	potential_new_points_and_following_data_directions.clear();
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> potential_new_points_data;
	bool direction1_valid, direction2_valid;
	compute_all_potential_3d_points_3views(sfmd, all_2d_reprojections_3views,all_2d_reprojections_ids_3views,potential_new_points_data);

	for(auto &potential_new_point : potential_new_points_data) {
		// check PLG compatibility
		vector<new_3dpoint_plgp_matches> valid_points_direction2,valid_points_direction1;
		vector<ulong> directions1, directions2;
		if(compatible_new_plg_point(sfmd, all_fundamental_matrices, plgs, potential_new_point, directions1, direction1_valid, valid_points_direction1, directions2, direction2_valid, valid_points_direction2)) {
			// 3D point has been considered compatible with PLG following
			potential_new_points_and_following_data_directions.push_back(make_new_3dpoint_and_sides_plgp_matches(valid_points_direction1,directions1,potential_new_point,valid_points_direction2,directions2));
		}
	}
}

void compute_unique_potential_3d_points_3views_plg_following_newpoint_compatibility(const SfMData &sfmd,const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs,
		const tuple<vector<PolyLineGraph2D::plg_point>,vector<PolyLineGraph2D::plg_point>,vector<PolyLineGraph2D::plg_point>> &all_2d_reprojections_3views,
		const tuple<int,int,int> &all_2d_reprojections_ids_3views,
		new_3dpoint_and_sides_plgp_matches &new_point, bool &valid) {
	valid = true;
	bool found = false;

	std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> potential_new_point;

	bool direction1_valid, direction2_valid;
	vector<new_3dpoint_plgp_matches> valid_points_direction2,valid_points_direction1;
	vector<ulong> directions1, directions2;

	bool tvalid;
	vector<int> views;
	views.push_back(get<0>(all_2d_reprojections_ids_3views));
	views.push_back(get<1>(all_2d_reprojections_ids_3views));
	views.push_back(get<2>(all_2d_reprojections_ids_3views));

	for(const auto &plgp0 : get<0>(all_2d_reprojections_3views))
		for(const auto &plgp1 : get<1>(all_2d_reprojections_3views))
			for(const auto &plgp2 : get<2>(all_2d_reprojections_3views)) {
				glm::vec3 triangulated_point;
				vector<PolyLineGraph2D::plg_point> new_plgps;
				new_plgps.push_back(plgp0);
				new_plgps.push_back(plgp1);
				new_plgps.push_back(plgp2);

				em_estimate3Dpositions(sfmd, new_plgps, views, triangulated_point, tvalid);

				vector<int> views_cpy = views;
				if(tvalid) {
					// Try following PLG from the new point

					potential_new_point=make_tuple(triangulated_point,new_plgps,views_cpy);

					if(compatible_new_plg_point(sfmd, all_fundamental_matrices, plgs, potential_new_point, directions1, direction1_valid, valid_points_direction1, directions2, direction2_valid, valid_points_direction2)) {
						if(found) {
							// Another point was already found, break
							valid = false;
							return;
						} else {
							found = true;
							new_point = make_new_3dpoint_and_sides_plgp_matches(valid_points_direction1,directions1,potential_new_point,valid_points_direction2,directions2);
						}
					}
				}
			}

	if(!found)
		valid=false;
}

void compute_findfirst_potential_3d_points_3views_plg_following_newpoint_compatibility(const SfMData &sfmd,const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs,
		const tuple<vector<PolyLineGraph2D::plg_point>,vector<PolyLineGraph2D::plg_point>,vector<PolyLineGraph2D::plg_point>> &all_2d_reprojections_3views,
		const tuple<int,int,int> &all_2d_reprojections_ids_3views,
		new_3dpoint_and_sides_plgp_matches &new_point, bool &valid) {
	valid = false;

	std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> potential_new_point;

	bool direction1_valid, direction2_valid;
	vector<new_3dpoint_plgp_matches> valid_points_direction2,valid_points_direction1;
	vector<ulong> directions1, directions2;

	bool tvalid;
	vector<int> views;
	views.push_back(get<0>(all_2d_reprojections_ids_3views));
	views.push_back(get<1>(all_2d_reprojections_ids_3views));
	views.push_back(get<2>(all_2d_reprojections_ids_3views));

	for(const auto &plgp0 : get<0>(all_2d_reprojections_3views))
		for(const auto &plgp1 : get<1>(all_2d_reprojections_3views))
			for(const auto &plgp2 : get<2>(all_2d_reprojections_3views)) {
				glm::vec3 triangulated_point;
				vector<PolyLineGraph2D::plg_point> new_plgps;
				new_plgps.push_back(plgp0);
				new_plgps.push_back(plgp1);
				new_plgps.push_back(plgp2);

				em_estimate3Dpositions(sfmd, new_plgps, views, triangulated_point, tvalid);

				vector<int> views_cpy = views;
				if(tvalid) {
					potential_new_point=make_tuple(triangulated_point,new_plgps,views_cpy);

					if(compatible_new_plg_point(sfmd, all_fundamental_matrices, plgs, potential_new_point, directions1, direction1_valid, valid_points_direction1, directions2, direction2_valid, valid_points_direction2)) {
						valid = true;
						new_point = make_new_3dpoint_and_sides_plgp_matches(valid_points_direction1,directions1,potential_new_point,valid_points_direction2,directions2);
						return;
					}
				}
			}
}

bool expand_point_to_other_view(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const int starting_plg_id, const int other_plg_id, const vector<PolyLineGraph2D::plg_point> &epipolar_correspondences, new_3dpoint_and_sides_plgp_matches &current_p3d_with_sides) {
	for(auto &epc : epipolar_correspondences)
		// try to add epc to current p3d
		if(add_view_to_3dpoint_and_sides_plgp_matches(sfmd, all_fundamental_matrices, plgs, current_p3d_with_sides, other_plg_id, epc))
			return true; // No adding more than one correspondence per image
	return false;
}

void expand_point_to_other_views(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const int starting_plg_id, const vector<vector<PolyLineGraph2D::plg_point>> &epipolar_correspondences, const int selected_views[], new_3dpoint_and_sides_plgp_matches &current_p3d_with_sides) {

	for(int i=0; i < selected_views[0]; i++)
		expand_point_to_other_view(sfmd, plgs, all_fundamental_matrices, starting_plg_id, i, epipolar_correspondences[i], current_p3d_with_sides);

	for(int i=selected_views[0]+1; i < selected_views[1]; i++)
		expand_point_to_other_view(sfmd, plgs, all_fundamental_matrices, starting_plg_id, i, epipolar_correspondences[i], current_p3d_with_sides);

	for(int i=selected_views[1]+1; i < selected_views[2]; i++)
		expand_point_to_other_view(sfmd, plgs, all_fundamental_matrices, starting_plg_id, i, epipolar_correspondences[i], current_p3d_with_sides);

	for(int i=selected_views[2]+1; i < sfmd.numCameras_; i++)
		expand_point_to_other_view(sfmd, plgs, all_fundamental_matrices, starting_plg_id, i, epipolar_correspondences[i], current_p3d_with_sides);

}

void compute_3D_point_multiple_views_plg_following(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const int starting_plg_id, const vector<vector<PolyLineGraph2D::plg_point>> &epipolar_correspondences, new_3dpoint_and_sides_plgp_matches &p3d_with_sides, bool &valid) {
	valid = false;

	int amount_of_non_empty = 0;
	int min_index,max_index;
	min_index = -1;
	max_index = -1;
	for(int i=0; i < epipolar_correspondences.size();i++)
		if(epipolar_correspondences[i].size()>0) {
			amount_of_non_empty++;
			min_index = min_index != -1 ? min_index : i;
			max_index = i;
		}

	if(amount_of_non_empty < 3)
		return;

	int cur_relative_index=0;
	int relative_mid_index = amount_of_non_empty/2;
	int mid_index;
	for(int i=0; i < epipolar_correspondences.size();i++)
		if(epipolar_correspondences[i].size()>0) {
			if(cur_relative_index == relative_mid_index) {
				mid_index = i;
				break;
			} else
				cur_relative_index++;
		}


	int selected_indexes[3] = {
			min_index,
			starting_plg_id == min_index || starting_plg_id == max_index ? mid_index : starting_plg_id,
			max_index
	};

	vector<new_3dpoint_and_sides_plgp_matches> p3ds_with_sides_3views;

	tuple<vector<PolyLineGraph2D::plg_point>,vector<PolyLineGraph2D::plg_point>,vector<PolyLineGraph2D::plg_point>> all_2d_reprojections_3views = make_tuple(epipolar_correspondences[selected_indexes[0]],epipolar_correspondences[selected_indexes[1]],epipolar_correspondences[selected_indexes[2]]);
	tuple<int,int,int> all_2d_reprojections_ids_3views = make_tuple(selected_indexes[0],selected_indexes[1],selected_indexes[2]);

	compute_unique_potential_3d_points_3views_plg_following_newpoint_compatibility(sfmd,all_fundamental_matrices, plgs,
				all_2d_reprojections_3views,
				all_2d_reprojections_ids_3views,
				p3d_with_sides,valid);

	if(!valid)
		return;

	expand_point_to_other_views(sfmd, plgs, all_fundamental_matrices, starting_plg_id, epipolar_correspondences, selected_indexes, p3d_with_sides);
}

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> compute_3D_point_multiple_views_plg_following_vecpoints(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const int starting_plg_id, const vector<vector<PolyLineGraph2D::plg_point>> &epipolar_correspondences) {
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> res;
	new_3dpoint_and_sides_plgp_matches p3d_with_sides;
	bool valid;

	compute_3D_point_multiple_views_plg_following(sfmd, plgs, all_fundamental_matrices, starting_plg_id, epipolar_correspondences, p3d_with_sides, valid);

	if(!valid)
		return res;

	for(auto &p3d : get<0>(p3d_with_sides).first)
		res.push_back(p3d);
	res.push_back(get<1>(p3d_with_sides));
	for(auto &p3d : get<2>(p3d_with_sides).first)
		res.push_back(p3d);

	return res;
}

#define SWITCH_DISABLE_INTERVAL
//#define TEST_SWITCH_EXPANDCENTRALPOINTONLY
void expand_allpoints_to_other_view_using_plmap(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const int starting_plg_id, const int other_plg_id, const vector<PolyLineGraph2D::plg_point> &epipolar_correspondences, const PolyLine2DMapSearch &plmap, vector<new_3dpoint_plgp_matches> &cur_p3ds, vector<ulong> &start_dirs, vector<ulong> &end_dirs, int &original_central_point_index) {
#if defined(SWITCH_PLG_MATCHING_ADDPOINT_BOTHDIR_ONE)

	bool success;

	pair<int,int> central_point_epc_matches;
	pair<int,int> central_point_epc_matches_indices;
	bool central_point_epc_matched=false;

	const int p3ds_initial_size = cur_p3ds.size();

	for(auto &epc : epipolar_correspondences) {
		// try to add epc to current p3d
		central_point_epc_matches = add_view_to_3dpoint_and_sides_plgp_matches_vector(sfmd, all_fundamental_matrices, plgs, cur_p3ds, start_dirs, end_dirs, other_plg_id, epc, 0, original_central_point_index, cur_p3ds.size(),central_point_epc_matched);
		if(central_point_epc_matched) {
			if(central_point_epc_matches.first > original_central_point_index) {
				// points added at start
				original_central_point_index = central_point_epc_matches.first;
				central_point_epc_matches_indices.first = 0;
				central_point_epc_matches_indices.second = central_point_epc_matches.first + central_point_epc_matches.second;
			} else {
				central_point_epc_matches_indices.first = original_central_point_index - central_point_epc_matches.first;
				central_point_epc_matches_indices.second = original_central_point_index + central_point_epc_matches.second;
			}
			break;
		}
	}

	int cur_interval_start=0;
	ulong cur_interval_pl_id;
	int start_match;

	ulong cur_pl_id;
	bool valid;

	pair<int,int> added_points_start_end;

	int next_point;
	glm::vec2 prev_coords;
	const glm::vec2 start_coords;

	int central_point;

	int cur_interval_end;

#if !defined(TEST_SWITCH_EXPANDCENTRALPOINTONLY)

#if defined(SWITCH_DISABLE_INTERVAL)
	int last_matched=-1;
	for(int cur_point=0; cur_point < cur_p3ds.size(); cur_point++) {
		if(central_point_epc_matched && cur_point == central_point_epc_matches_indices.first) {
			cur_point = central_point_epc_matches_indices.second; // Skip to after the epc interval
			last_matched = central_point_epc_matches_indices.second;
			continue;
		}

		start_coords = compute_projection(sfmd,other_plg_id,get<0>(cur_p3ds[cur_point]));
		plmap.find_unique_polyline_potentially_within_search_dist(start_coords,cur_pl_id,valid);

		if(valid) {
			central_point = cur_point;

			// Compute PLGP on new view, polyline cur_interval_pl_id
			const PolyLineGraph2D::polyline &pl = plgs[other_plg_id].polylines[cur_pl_id];
			PolyLineGraph2D::polyline::pl_point init_ppl;
			if(pl.compute_distancesq(start_coords,init_ppl) > MAX_3DPOINT_PROJECTIONDISTSQ_EXPANDALLVIEWS)
				return false;
			PolyLineGraph2D::plg_point init_plgp(cur_pl_id,init_ppl);

			if(central_point_epc_matched)
				cur_interval_end = central_point <= central_point_epc_matches_indices.first ? central_point_epc_matches_indices.first : cur_p3ds.size();
			else
				cur_interval_end = cur_p3ds.size();

			// Try to add this point, get last matched index
			added_points_start_end = add_view_to_3dpoint_and_sides_plgp_matches_vector(sfmd, all_fundamental_matrices, plgs, cur_p3ds, start_dirs, end_dirs, other_plg_id, init_plgp, last_matched+1, central_point, cur_interval_end,success);

			if(success) {
				if(added_points_start_end.first > central_point) {
					original_central_point_index = added_points_start_end.first;
					// Points were added before start
					cur_point = (added_points_start_end.first+ added_points_start_end.second);
				} else {
					cur_point = central_point + (added_points_start_end.second);
				}
				last_matched = cur_point;
			}	
		}
	}
#else
	for(int cur_point=0; cur_point < cur_p3ds.size(); cur_point++) {
		if(central_point_epc_matched && cur_point == central_point_epc_matches_indices.first) {
			cur_point = central_point_epc_matches_indices.second; // Skip to after the epc interval

			// Reset interval
			cur_interval_start=cur_point+1;
			continue;
		}

		start_coords = compute_projection(sfmd,other_plg_id,get<0>(cur_p3ds[cur_point]));
		plmap.find_unique_polyline_potentially_within_search_dist(start_coords,cur_pl_id,valid);

		if(!valid || (cur_interval_start != cur_point && cur_interval_pl_id != cur_pl_id)) {
			// Restart interval
			cur_interval_start=cur_point+1;
			//cur_interval_pl_id=cur_pl_id;
		}
		else
		{
			if(cur_interval_start == cur_point) {
				cur_interval_pl_id=cur_pl_id;
			} else if(cur_point - cur_interval_start == 1) {
			// A sequence of 3 valid points has been found => it's time to run plg matching

			central_point = cur_point;

			// Compute PLGP on new view, polyline cur_interval_pl_id
			const PolyLineGraph2D::polyline &pl = plgs[other_plg_id].polylines[cur_interval_pl_id];
			PolyLineGraph2D::polyline::pl_point init_ppl;
			if(pl.compute_distancesq(prev_coords,init_ppl) > MAX_3DPOINT_PROJECTIONDISTSQ_EXPANDALLVIEWS)
				return false;
			PolyLineGraph2D::plg_point init_plgp(cur_interval_pl_id,init_ppl);

			if(central_point_epc_matched)
				cur_interval_end = central_point <= central_point_epc_matches_indices.first ? central_point_epc_matches_indices.first : cur_p3ds.size();
			else
				cur_interval_end = cur_p3ds.size();

			// Try to add this point, get last matched index
			added_points_start_end = add_view_to_3dpoint_and_sides_plgp_matches_vector(sfmd, all_fundamental_matrices, plgs, cur_p3ds, start_dirs, end_dirs, other_plg_id, init_plgp, cur_interval_start, central_point, cur_interval_end,success);

			if(added_points_start_end.first > central_point) {
				original_central_point_index = added_points_start_end.first;
				// Points were added before start
				cur_point = (added_points_start_end.first+ added_points_start_end.second);
			} else
				cur_point = central_point + (added_points_start_end.second);

			// Reset interval
			cur_interval_start = cur_point + 1;
		}

		prev_coords = start_coords;
	}
	}
#endif
#else
	if(!central_point_epc_matched) {
		int cur_point = original_central_point_index;

		start_coords = compute_projection(sfmd,other_plg_id,get<0>(cur_p3ds[cur_point]));
		plmap.find_unique_polyline_potentially_within_search_dist(start_coords,cur_pl_id,valid);

		if(valid) {
			// A sequence of 3 valid points has been found => it's time to run plg matching

			central_point = cur_point;

			// Compute PLGP on new view, polyline cur_interval_pl_id
			const PolyLineGraph2D::polyline &pl = plgs[other_plg_id].polylines[cur_pl_id];
			PolyLineGraph2D::polyline::pl_point init_ppl;
			if(pl.compute_distancesq(start_coords,init_ppl) > MAX_3DPOINT_PROJECTIONDISTSQ_EXPANDALLVIEWS)
				return false;
			PolyLineGraph2D::plg_point init_plgp(cur_pl_id,init_ppl);

			cur_interval_end = cur_p3ds.size();

			// Try to add this point, get last matched index
			added_points_start_end = add_view_to_3dpoint_and_sides_plgp_matches_vector(sfmd, all_fundamental_matrices, plgs, cur_p3ds, start_dirs, end_dirs, other_plg_id, init_plgp, 0, central_point, cur_interval_end,success);

			if(added_points_start_end.first > central_point)
				original_central_point_index = added_points_start_end.first;
		}
	}
#endif

#endif
}

bool expand_point_to_other_view_using_plmap(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const int starting_plg_id, const int other_plg_id, const PolyLine2DMapSearch &plmap, new_3dpoint_and_sides_plgp_matches &current_p3d_with_sides) {
	const new_3dpoint_plgp_matches &central_p3d = get<1>(current_p3d_with_sides);
	const glm::vec2 &start_coords = compute_projection(sfmd,other_plg_id,get<0>(central_p3d));
	ulong pl_id;
	bool valid;
	plmap.find_unique_polyline_potentially_within_search_dist(start_coords,pl_id,valid);
	if(valid) {
		const PolyLineGraph2D::polyline &pl = plgs[other_plg_id].polylines[pl_id];
		PolyLineGraph2D::polyline::pl_point init_ppl;
		if(pl.compute_distancesq(start_coords,init_ppl) > MAX_3DPOINT_PROJECTIONDISTSQ_EXPANDALLVIEWS)
			return false;
		PolyLineGraph2D::plg_point init_plgp(pl_id,init_ppl);
		if(add_view_to_3dpoint_and_sides_plgp_matches(sfmd, all_fundamental_matrices, plgs, current_p3d_with_sides, other_plg_id, init_plgp)) {
			//cout << "Successfully added one more view through plmap" << endl;
			return true; // No adding more than one correspondence per image
		}
	}
	return false;
}

void expand_point_to_other_views_expandallviews(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const int starting_plg_id, const vector<vector<PolyLineGraph2D::plg_point>> &epipolar_correspondences, const int selected_views[], const vector<PolyLine2DMapSearch> &plmaps, new_3dpoint_and_sides_plgp_matches &current_p3d_with_sides) {

	for(int i=0; i < selected_views[0]; i++)
		if(!expand_point_to_other_view(sfmd, plgs, all_fundamental_matrices, starting_plg_id, i, epipolar_correspondences[i], current_p3d_with_sides))
			expand_point_to_other_view_using_plmap(sfmd, plgs, all_fundamental_matrices, starting_plg_id, i, plmaps[i], current_p3d_with_sides);

	for(int i=selected_views[0]+1; i < selected_views[1]; i++)
		if(!expand_point_to_other_view(sfmd, plgs, all_fundamental_matrices, starting_plg_id, i, epipolar_correspondences[i], current_p3d_with_sides))
			expand_point_to_other_view_using_plmap(sfmd, plgs, all_fundamental_matrices, starting_plg_id, i, plmaps[i], current_p3d_with_sides);

	for(int i=selected_views[1]+1; i < selected_views[2]; i++)
		if(!expand_point_to_other_view(sfmd, plgs, all_fundamental_matrices, starting_plg_id, i, epipolar_correspondences[i], current_p3d_with_sides))
			expand_point_to_other_view_using_plmap(sfmd, plgs, all_fundamental_matrices, starting_plg_id, i, plmaps[i], current_p3d_with_sides);

	for(int i=selected_views[2]+1; i < sfmd.numCameras_; i++)
		if(!expand_point_to_other_view(sfmd, plgs, all_fundamental_matrices, starting_plg_id, i, epipolar_correspondences[i], current_p3d_with_sides))
			expand_point_to_other_view_using_plmap(sfmd, plgs, all_fundamental_matrices, starting_plg_id, i, plmaps[i], current_p3d_with_sides);
}

void expand_point_to_other_views_expandallviews_vector(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const int starting_plg_id, const vector<vector<PolyLineGraph2D::plg_point>> &epipolar_correspondences, const int selected_views[], const vector<PolyLine2DMapSearch> &plmaps, vector<new_3dpoint_plgp_matches> &cur_p3ds, vector<ulong> &start_dirs, vector<ulong> &end_dirs, int &central_point_index) {

	for(int i=0; i < selected_views[0]; i++)
		expand_allpoints_to_other_view_using_plmap(sfmd, plgs, all_fundamental_matrices, starting_plg_id, i, epipolar_correspondences[i], plmaps[i], cur_p3ds, start_dirs, end_dirs, central_point_index);

	for(int i=selected_views[0]+1; i < selected_views[1]; i++)
		expand_allpoints_to_other_view_using_plmap(sfmd, plgs, all_fundamental_matrices, starting_plg_id, i, epipolar_correspondences[i], plmaps[i], cur_p3ds, start_dirs, end_dirs, central_point_index);

	for(int i=selected_views[1]+1; i < selected_views[2]; i++)
		expand_allpoints_to_other_view_using_plmap(sfmd, plgs, all_fundamental_matrices, starting_plg_id, i, epipolar_correspondences[i], plmaps[i], cur_p3ds, start_dirs, end_dirs, central_point_index);

	for(int i=selected_views[2]+1; i < sfmd.numCameras_; i++)
		expand_allpoints_to_other_view_using_plmap(sfmd, plgs, all_fundamental_matrices, starting_plg_id, i, epipolar_correspondences[i], plmaps[i], cur_p3ds, start_dirs, end_dirs, central_point_index);
}

void compute_3D_point_multiple_views_plg_following_expandallviews(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const int starting_plg_id, const vector<vector<PolyLineGraph2D::plg_point>> &epipolar_correspondences, const vector<PolyLine2DMapSearch> &plmaps, new_3dpoint_and_sides_plgp_matches &p3d_with_sides, bool &valid) {
	valid = false;

	int amount_of_non_empty = 0;
	int min_index,max_index;
	min_index = -1;
	max_index = -1;
	for(int i=0; i < epipolar_correspondences.size();i++)
		if(epipolar_correspondences[i].size()>0) {
			amount_of_non_empty++;
			min_index = min_index != -1 ? min_index : i;
			max_index = i;
		}

	if(amount_of_non_empty < 3)
		return;

	int cur_relative_index=0;
	int relative_mid_index = amount_of_non_empty/2;
	int mid_index;
	for(int i=0; i < epipolar_correspondences.size();i++)
		if(epipolar_correspondences[i].size()>0) {
			if(cur_relative_index == relative_mid_index) {
				mid_index = i;
				break;
			} else
				cur_relative_index++;
		}


	int selected_indexes[3] = {
			min_index,
			starting_plg_id == min_index || starting_plg_id == max_index ? mid_index : starting_plg_id,
			max_index
	};

	vector<new_3dpoint_and_sides_plgp_matches> p3ds_with_sides_3views;

	tuple<vector<PolyLineGraph2D::plg_point>,vector<PolyLineGraph2D::plg_point>,vector<PolyLineGraph2D::plg_point>> all_2d_reprojections_3views = make_tuple(epipolar_correspondences[selected_indexes[0]],epipolar_correspondences[selected_indexes[1]],epipolar_correspondences[selected_indexes[2]]);
	tuple<int,int,int> all_2d_reprojections_ids_3views = make_tuple(selected_indexes[0],selected_indexes[1],selected_indexes[2]);

	compute_unique_potential_3d_points_3views_plg_following_newpoint_compatibility(sfmd,all_fundamental_matrices, plgs,
				all_2d_reprojections_3views,
				all_2d_reprojections_ids_3views,
				p3d_with_sides,valid);

	if(!valid)
		return;

	expand_point_to_other_views_expandallviews(sfmd, plgs, all_fundamental_matrices, starting_plg_id, epipolar_correspondences, selected_indexes, plmaps, p3d_with_sides);
}

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> compute_3D_point_multiple_views_plg_following_expandallviews_vector(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const int starting_plg_id, const vector<vector<PolyLineGraph2D::plg_point>> &epipolar_correspondences, const vector<PolyLine2DMapSearch> &plmaps) {
	vector<new_3dpoint_plgp_matches> res;

	new_3dpoint_and_sides_plgp_matches p3d_with_sides;
	bool valid;

	valid = false;

	int amount_of_non_empty = 0;
	int min_index,max_index;
	min_index = -1;
	max_index = -1;
	for(int i=0; i < epipolar_correspondences.size();i++)
		if(epipolar_correspondences[i].size()>0) {
			amount_of_non_empty++;
			min_index = min_index != -1 ? min_index : i;
			max_index = i;
		}

	if(amount_of_non_empty < 3)
		return res;

	int cur_relative_index=0;
	int relative_mid_index = amount_of_non_empty/2;
	int mid_index;
	for(int i=0; i < epipolar_correspondences.size();i++)
		if(epipolar_correspondences[i].size()>0) {
			if(cur_relative_index == relative_mid_index) {
				mid_index = i;
				break;
			} else
				cur_relative_index++;
		}


	int selected_indexes[3] = {
			min_index,
			starting_plg_id == min_index || starting_plg_id == max_index ? mid_index : starting_plg_id,
			max_index
	};

	vector<new_3dpoint_and_sides_plgp_matches> p3ds_with_sides_3views;

	tuple<vector<PolyLineGraph2D::plg_point>,vector<PolyLineGraph2D::plg_point>,vector<PolyLineGraph2D::plg_point>> all_2d_reprojections_3views = make_tuple(epipolar_correspondences[selected_indexes[0]],epipolar_correspondences[selected_indexes[1]],epipolar_correspondences[selected_indexes[2]]);
	tuple<int,int,int> all_2d_reprojections_ids_3views = make_tuple(selected_indexes[0],selected_indexes[1],selected_indexes[2]);

	compute_unique_potential_3d_points_3views_plg_following_newpoint_compatibility(sfmd,all_fundamental_matrices, plgs,
				all_2d_reprojections_3views,
				all_2d_reprojections_ids_3views,
				p3d_with_sides,valid);

	if(!valid)
		return res;

	res = new_3dpoint_and_sides_plgp_matches_to_vector(p3d_with_sides);

	int central_point = get<0>(p3d_with_sides).first.size();

	expand_point_to_other_views_expandallviews_vector(sfmd, plgs, all_fundamental_matrices, starting_plg_id, epipolar_correspondences, selected_indexes, plmaps, res, get<0>(p3d_with_sides).second, get<2>(p3d_with_sides).second, central_point);

	return res;
}

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> compute_3D_point_multiple_views_plg_following_vecpoints_expandallviews(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const int starting_plg_id, const vector<vector<PolyLineGraph2D::plg_point>> &epipolar_correspondences, const vector<PolyLine2DMapSearch> &plmaps) {
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> res;
	new_3dpoint_and_sides_plgp_matches p3d_with_sides;
	bool valid;

	compute_3D_point_multiple_views_plg_following_expandallviews(sfmd, plgs, all_fundamental_matrices, starting_plg_id, epipolar_correspondences, plmaps, p3d_with_sides, valid);

	if(!valid)
		return res;

	res = new_3dpoint_and_sides_plgp_matches_to_vector(p3d_with_sides);

	return res;
}

void compute_3d_point_coords_combinations(const SfMData &sfmd,
		const vector<glm::vec2> &all_2d_reprojections_coords,
		const vector<int> &all_2d_reprojections_ids,
		const int min_combinations,
		vector<glm::vec2> &selected_2d_reprojections_coords,
		vector<int> &selected_2d_reprojections_ids,
		vector<bool> &selected,
		glm::vec3 &new_point_data,
		bool &valid) {
	valid = false;

	selected.resize(all_2d_reprojections_ids.size());
	std::fill(selected.begin() + min_combinations, selected.end(), false);
    std::fill(selected.begin(), selected.begin() + min_combinations, true);

    do {
    	selected_2d_reprojections_coords.clear();
    	selected_2d_reprojections_ids.clear();

        for (int i = 0; i < all_2d_reprojections_ids.size(); ++i) {
            if (selected[i]) {
            	selected_2d_reprojections_coords.push_back(all_2d_reprojections_coords[i]);
            	selected_2d_reprojections_ids.push_back(all_2d_reprojections_ids[i]);
            }
        }

        compute_3d_point_coords(sfmd,selected_2d_reprojections_coords,selected_2d_reprojections_ids,new_point_data,valid);
    } while (!valid && std::prev_permutation(selected.begin(), selected.end()));

    if(!valid)
    	return;

    glm::vec3 new_3d_coords;
    bool need_reorder=false;
    for (int i = 0; i < all_2d_reprojections_ids.size(); ++i) {
        if (!selected[i]) {
        	em_add_new_observation_to_3Dpositions(sfmd, new_point_data, selected_2d_reprojections_coords, selected_2d_reprojections_ids, all_2d_reprojections_coords[i], all_2d_reprojections_ids[i], new_3d_coords, valid);
        	if(valid) {
        		selected[i] = true;
        		new_point_data = new_3d_coords;
        		selected_2d_reprojections_coords.push_back(all_2d_reprojections_coords[i]);
        		if(!need_reorder && all_2d_reprojections_ids[i] > selected_2d_reprojections_ids[selected_2d_reprojections_ids.size()-1])
        			need_reorder = true;
        		selected_2d_reprojections_ids.push_back(all_2d_reprojections_ids[i]);
        	}
        }
    }

    valid = true;

    // Re-Order
    if(need_reorder)
    	reorder_pair_of_vector(selected_2d_reprojections_coords,selected_2d_reprojections_ids);
}

void compute_3d_point(const SfMData &sfmd,
		const vector<glm::vec2> &selected_2d_reprojections_coords,
		const vector<int> &selected_2d_reprojections_ids,
		std::tuple<glm::vec3, vector<glm::vec2>, vector<int>> &new_point_data,
		bool &valid) {
	valid = false;

	if (selected_2d_reprojections_coords.size() >= 2) {
		em_estimate3Dpositions(sfmd, selected_2d_reprojections_coords,
				selected_2d_reprojections_ids, std::get < 0 > (new_point_data),
				valid);
		if (valid) {
			std::get < 1 > (new_point_data) = selected_2d_reprojections_coords;
			std::get < 2 > (new_point_data) = selected_2d_reprojections_ids;
		}
	}
}
