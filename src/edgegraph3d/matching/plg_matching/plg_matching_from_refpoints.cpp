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


#include "plg_matching_from_refpoints.hpp"

#include <omp.h>
#include <opencv2/core/mat.hpp>
#include <iostream>
#include <utility>
#include <vector>

#include "polyline_graph_2d_hmap_impl.hpp"
#include "SfMData.h"
#include "plg_edge_manager.hpp"
#include "plg_edge_manager_closest_only.hpp"
#include "plgp_consensus_manager.hpp"
#include "plg_matches_manager.hpp"
#include "plg_matching.hpp"
#include "polyline_graph_2d.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "global_defines.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"
#include "glm.hpp"

#include "drawing_utilities.hpp"

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> plg_matching_from_refpoint_starting_image(const SfMData &sfm_data, const EdgeManager *em,const PLGPConsensusManager *cm, const ulong refpoint_id, const ulong starting_img_id) {
	std::pair<std::vector<PolyLineGraph2D::plg_point>,std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>>> intersections_and_correspondences_pair = ((PLGEdgeManagerClosestOnly*) em)->detect_nearby_intersections_and_correspondences_plgp(starting_img_id, refpoint_id, DETECTION_STARTING_RADIUS, DETECTION_CORRESPONDENCES_MULTIPLICATION_FACTOR);

	std::vector<PolyLineGraph2D::plg_point> &intersections = intersections_and_correspondences_pair.first;
	std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>> &correspondences = intersections_and_correspondences_pair.second;

	return cm->consensus_strategy_single_point(starting_img_id, refpoint_id, intersections_and_correspondences_pair);
}

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> plg_matching_from_refpoint(const SfMData &sfm_data, const EdgeManager *em,const PLGPConsensusManager *cm, const ulong refpoint_id, PLGMatchesManager &plgmm) {
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> res;

	std::vector<std::pair<std::vector<PolyLineGraph2D::plg_point>,std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>>>> intersections_and_correspondences_all_imgs = ((PLGEdgeManager*) em)->detect_nearby_intersections_and_correspondences_plgp(refpoint_id);

	for(int i=0; i < sfm_data.camViewingPointN_[refpoint_id].size(); i++) {
		int starting_img_id = sfm_data.camViewingPointN_[refpoint_id][i];

		vector<vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>>> cur_res_vec = cm->consensus_strategy_single_point_vector(starting_img_id, refpoint_id, intersections_and_correspondences_all_imgs[i]);

		for(const auto &cur_res: cur_res_vec) {
			plgmm.add_matched_3dpolyline(cur_res);
			res.insert(res.end(),cur_res.begin(),cur_res.end());
		}
	}

	return res;
}

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> plg_matching_from_refpoints_parallel(const SfMData &sfm_data, const EdgeManager *em,const PLGPConsensusManager *cm, PLGMatchesManager &plgmm) {
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> res;

	int num_threads = omp_get_max_threads();
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> *res_t = new vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>>[num_threads];

#pragma omp for
	for(ulong refpoint_id=0; refpoint_id < sfm_data.numPoints_; refpoint_id++) {
		int thread_id = omp_get_thread_num();
		cout << "Extracting 3D edges from refpoint " << refpoint_id << "\n";
		vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> cur_res = plg_matching_from_refpoint(sfm_data, em, cm, refpoint_id,plgmm);
		res_t[thread_id].insert(res_t[thread_id].end(),cur_res.begin(),cur_res.end());
	}

	for(int i=0; i < num_threads; i++)
		for(auto &new_p3d : res_t[i])
			res.push_back(new_p3d);

	delete[] res_t;

	return res;
}

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> plg_matching_from_refpoints(const SfMData &sfm_data, const EdgeManager *em,const PLGPConsensusManager *cm, PLGMatchesManager &plgmm) {
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> res;

	for(ulong refpoint_id=0; refpoint_id < sfm_data.numPoints_; refpoint_id++) {
		cout << "Extracting 3D edges from refpoint " << refpoint_id << "\n";
		vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> cur_res = plg_matching_from_refpoint(sfm_data, em, cm, refpoint_id,plgmm);
		res.insert(res.end(),cur_res.begin(),cur_res.end());
	}

	return res;
}

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> plg_matching_from_refpoint(const SfMData &sfm_data, const EdgeManager *em,const PLGPConsensusManager *cm, const ulong refpoint_id) {
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> res;

	std::vector<std::pair<std::vector<PolyLineGraph2D::plg_point>,std::vector<std::vector<std::vector<PolyLineGraph2D::plg_point>>>>> intersections_and_correspondences_all_imgs = ((PLGEdgeManager*) em)->detect_nearby_intersections_and_correspondences_plgp(refpoint_id);

	for(int i=0; i < sfm_data.camViewingPointN_[refpoint_id].size(); i++) {
		int starting_img_id = sfm_data.camViewingPointN_[refpoint_id][i];
		vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> cur_res = cm->consensus_strategy_single_point(starting_img_id, refpoint_id, intersections_and_correspondences_all_imgs[i]);
		res.insert(res.end(),cur_res.begin(),cur_res.end());
	}

	return res;
}

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> plg_matching_from_refpoints_parallel(const SfMData &sfm_data, const EdgeManager *em,const PLGPConsensusManager *cm) {
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> res;

	int num_threads = omp_get_max_threads();
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> *res_t = new vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>>[num_threads];

#pragma omp for
	for(ulong refpoint_id=0; refpoint_id < sfm_data.numPoints_; refpoint_id++) {
		int thread_id = omp_get_thread_num();
		cout << "Extracting 3D edges from refpoint " << refpoint_id << "\n";
		vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> cur_res = plg_matching_from_refpoint(sfm_data, em, cm, refpoint_id);
		res_t[thread_id].insert(res_t[thread_id].end(),cur_res.begin(),cur_res.end());
	}

	for(int i=0; i < num_threads; i++)
		for(auto &new_p3d : res_t[i])
			res.push_back(new_p3d);

	delete[] res_t;

	return res;
}

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> plg_matching_from_refpoints(const SfMData &sfm_data, const EdgeManager *em,const PLGPConsensusManager *cm) {
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> res;

	for(ulong refpoint_id=0; refpoint_id < sfm_data.numPoints_; refpoint_id++) {
		cout << "Extracting 3D edges from refpoint " << refpoint_id << "\n";
		vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> cur_res = plg_matching_from_refpoint(sfm_data, em, cm, refpoint_id);
		res.insert(res.end(),cur_res.begin(),cur_res.end());
	}

	return res;
}

//pair<vector<std::tuple<glm::vec3,vector<glm::vec2>,vector<int>>>,tuple<vector<glm::vec3,PolyLineGraph2D::plg_point>,vector<int>>> process_refpoint(const vector<Mat> &imgs, vector<Mat> &out_imgs, const vector<PolyLineGraph2DHMapImpl> &plgs, const SfMData &sfmd, const Mat** all_fundamental_matrices, const PLGEdgeManager *em, const ulong refpoint_id, char * outfolder)
vector<std::tuple<glm::vec3,vector<glm::vec2>,vector<int>>> process_refpoint(const vector<Mat> &imgs, vector<Mat> &out_imgs, const vector<PolyLineGraph2DHMapImpl> &plgs, const SfMData &sfmd, const Mat** all_fundamental_matrices, const PLGEdgeManagerClosestOnly *em, const ulong refpoint_id, char * outfolder)
{
	cout << "Extracting 3D edges from refpoint " << refpoint_id << endl;

	vector<std::tuple<glm::vec3,vector<glm::vec2>,vector<int>>> res;
	const vector<int> &cams = sfmd.camViewingPointN_[refpoint_id];

	int starting_img_id;
	int starting_img_id_index;
	std::vector<pair<PolyLineGraph2D::plg_point,std::vector<std::vector<PolyLineGraph2D::plg_point>>>> nearby_potential_matches;// = ((PLGEdgeManager*) em)->detect_nearby_intersections_and_correspondences_plgp_exclude_parallel_epipsegments(starting_img_id, refpoint_id, DETECTION_STARTING_RADIUS, DETECTION_CORRESPONDENCES_MULTIPLICATION_FACTOR, DETECTION_CORRESPONDENCES_MIN_RADIUS, DETECTION_CORRESPONDENCES_DETECTION_MAX_ANGLECOS);

	for(starting_img_id_index = 0; starting_img_id_index < cams.size(); starting_img_id_index++) {
		starting_img_id = cams[starting_img_id_index];

		nearby_potential_matches = ((PLGEdgeManagerClosestOnly*) em)->detect_nearby_intersections_and_correspondences_plgp_exclude_parallel_epipsegments(starting_img_id, refpoint_id, DETECTION_STARTING_RADIUS, DETECTION_CORRESPONDENCES_MULTIPLICATION_FACTOR, DETECTION_CORRESPONDENCES_MIN_RADIUS, DETECTION_CORRESPONDENCES_DETECTION_MAX_ANGLECOS);

		if(nearby_potential_matches.size() > 0) {
			if(one_or_zero_correspondences(nearby_potential_matches[0]) && amount_of_total_2d_correspondences(nearby_potential_matches[0]) >= 3)
				break;
		}
	}

	if(starting_img_id_index >= cams.size())
		return res;

	const pair<PolyLineGraph2D::plg_point,std::vector<std::vector<PolyLineGraph2D::plg_point>>> &potential_match = nearby_potential_matches[0];

	pair<vector<glm::vec2>,vector<int>> p = convert_plgpoint_correspondences(cams, potential_match.second);
	std::tuple<glm::vec3, vector<glm::vec2>, vector<int>> new_point;
	bool valid;

	compute_3d_point(sfmd,p.first,p.second,new_point,valid);

	if(!valid)
		return res;

	pair<vector<PolyLineGraph2D::plg_point>,vector<int>> p_plgp = convert_plgpoint_correspondences_plgp(cams, potential_match.second);
	std::tuple<glm::vec3,vector<PolyLineGraph2D::plg_point>,vector<int>> new_point_plgp = make_tuple(get<0>(new_point),p_plgp.first,p.second);

	res.push_back(new_point);

	pair<vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>>,vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>>> new_points = follow_plgs_from_match3(sfmd, all_fundamental_matrices, plgs, new_point_plgp, valid);
	for(const auto &fp : new_points.first)
		res.push_back(make_tuple(get<0>(fp),convert_vecplgp_to_vec2(get<1>(fp)),get<2>(fp)));
	for(const auto &fp : new_points.second)
		res.push_back(make_tuple(get<0>(fp),convert_vecplgp_to_vec2(get<1>(fp)),get<2>(fp)));

    vector<Mat> cpy_imgs = copy_imgs(imgs);
	draw_3dpoints_on_imgs(cpy_imgs, res);
	draw_single_point_process(cpy_imgs,sfmd, all_fundamental_matrices, starting_img_id, refpoint_id, DETECTION_STARTING_RADIUS, nearby_potential_matches);

	draw_single_point_process_no_epilines(out_imgs,sfmd, all_fundamental_matrices, starting_img_id, refpoint_id, DETECTION_STARTING_RADIUS, nearby_potential_matches);

	return res;
}
