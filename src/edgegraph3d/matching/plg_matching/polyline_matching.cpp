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


#include "polyline_matching.hpp"

#include <omp.h>
#include <opencv2/core/mat.hpp>
#include <map>

#include "SfMData.h"
#include "plg_matches_manager.hpp"
#include "geometric_utilities.hpp"
#include "triangulation.hpp"
#include "global_switches.hpp"
#include "glm.hpp"


vector<vector<PolyLineGraph2D::plg_point>> find_epipolar_correspondences(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const vector<set<ulong>> &potentially_compatible_polylines, const int starting_plg_id, const PolyLineGraph2D::plg_point &starting_plgp, PLGMatchesManager &plgmm)
{
	PolyLineGraph2D::polyline::pl_interval pli_containing_point;
	vector<vector<PolyLineGraph2D::plg_point>> res;
	glm::vec3 epipolar;

	for(int other_plg_id=0; other_plg_id < plgs.size(); other_plg_id++) {
		vector<PolyLineGraph2D::plg_point> cur_res;
		vector<PolyLineGraph2D::plg_point> filtered_res;
		if(other_plg_id == starting_plg_id) {
			filtered_res.push_back(starting_plgp);
		} else {
			if(computeCorrespondEpilineSinglePoint(starting_plgp.plp.coords, all_fundamental_matrices[starting_plg_id][other_plg_id], epipolar, 1))
				for(set<ulong>::iterator jt=potentially_compatible_polylines[other_plg_id].begin(); jt != potentially_compatible_polylines[other_plg_id].end(); jt++)
				{
					ulong other_polyline_id = *jt;
					vector<PolyLineGraph2D::polyline::pl_point> polyline_points = plgs[other_plg_id].polylines[other_polyline_id].intersect_line(epipolar);
					cur_res = convert_vec_pl_point_to_plg_point(polyline_points,other_polyline_id);

					for(const auto &plgp : cur_res)
						if(!plgmm.is_matched(other_plg_id,plgp,pli_containing_point))
							filtered_res.push_back(plgp);
				}
		}
		res.push_back(filtered_res);
	}

	return res;
}



vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> find_new_3d_points_from_compatible_polylines_starting_plgp(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const vector<set<ulong>> &potentially_compatible_polylines, const int starting_plg_id, const PolyLineGraph2D::plg_point &starting_plgp, PLGMatchesManager &plgmm) {


	vector<vector<PolyLineGraph2D::plg_point>> epipolar_correspondences = find_epipolar_correspondences(sfmd, plgs, all_fundamental_matrices, potentially_compatible_polylines, starting_plg_id, starting_plgp,plgmm);

	return compute_3D_point_multiple_views_plg_following_vecpoints(sfmd, plgs, all_fundamental_matrices, starting_plg_id, epipolar_correspondences);
}


/**
 * Given
 * vector<vector<ulong>> : for each image, vector of polyline ids
 * representing sets of polylines on all images considered potentially compatible among them
 *
 * returns new 3D points vector
 */
vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> find_new_3d_points_from_compatible_polylines(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const vector<set<ulong>> &potentially_compatible_polylines, PLGMatchesManager &plgmm)
{
	vector<PolyLine2DMapSearch> plmaps;

	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> res;
	PolyLineGraph2D::polyline::pl_interval pli_containing_point;

	for(int starting_plg_id=0; starting_plg_id < plgs.size(); starting_plg_id++) {
		//cout << "Processing image " << starting_plg_id << "\n";
		for(const auto starting_polyline_id : potentially_compatible_polylines[starting_plg_id])
		{
			//cout << "Processing polyline " << starting_polyline_id << "\n";
			const PolyLineGraph2D::polyline &pl = plgs[starting_plg_id].polylines[starting_polyline_id];
			PolyLineGraph2D::polyline::pl_point plp = pl.get_start_plp();
			bool reached_end;
			plp = pl.next_pl_point_by_distance(plp,pl.end,SPLIT_INTERVAL_DISTANCE,reached_end);
			while(!reached_end) {
				//cout << "PLP start: " << plp.segment_index << " - " << plp.coords << "\n";
				if(!plgmm.is_matched(starting_plg_id,starting_polyline_id,plp,pli_containing_point)) {
					PolyLineGraph2D::plg_point starting_plgp(starting_polyline_id,plp);

					// Detect new 3D points using starting_plgp and potentially compatible polylines
					vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> cur_res = find_new_3d_points_from_compatible_polylines_starting_plgp(sfmd, plgs, all_fundamental_matrices, potentially_compatible_polylines, starting_plg_id, starting_plgp,plgmm);

					plgmm.add_matched_3dpolyline(cur_res);

					for(auto &new_p3d : cur_res)
						res.push_back(new_p3d);

					plp = pl.next_pl_point_by_distance(plp,pl.end,SPLIT_INTERVAL_DISTANCE,reached_end);
				} else {
					// jump to interval end
					plp = pl.next_pl_point_by_distance(pli_containing_point.end,pl.end,SPLIT_INTERVAL_DISTANCE,reached_end);
				}
			}
		}
	}

	return res;
}

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> find_new_3d_points_from_compatible_polylines_starting_plgp_expandallviews(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const vector<set<ulong>> &potentially_compatible_polylines, const int starting_plg_id, const PolyLineGraph2D::plg_point &starting_plgp, PLGMatchesManager &plgmm, const vector<PolyLine2DMapSearch> &plmaps) {


	vector<vector<PolyLineGraph2D::plg_point>> epipolar_correspondences = find_epipolar_correspondences(sfmd, plgs, all_fundamental_matrices, potentially_compatible_polylines, starting_plg_id, starting_plgp,plgmm);

#if defined(GLOBAL_SWITCH_USEEXPANDALLVIEWSVECTOR)
	return compute_3D_point_multiple_views_plg_following_expandallviews_vector(sfmd, plgs, all_fundamental_matrices, starting_plg_id, epipolar_correspondences,plmaps);
#else
	return compute_3D_point_multiple_views_plg_following_vecpoints_expandallviews(sfmd, plgs, all_fundamental_matrices, starting_plg_id, epipolar_correspondences,plmaps);
#endif
}

/**
 * Given
 * vector<vector<ulong>> : for each image, vector of polyline ids
 * representing sets of polylines on all images considered potentially compatible among them
 *
 * returns new 3D points vector
 */
vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> find_new_3d_points_from_compatible_polylines_expandallviews_parallel(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const vector<set<ulong>> &potentially_compatible_polylines, PLGMatchesManager &plgmm, const vector<PolyLine2DMapSearch> &plmaps)
{
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> res;
	PolyLineGraph2D::polyline::pl_interval pli_containing_point;

	int num_threads = omp_get_max_threads();
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> *res_t = new vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>>[num_threads];

#pragma omp for private(pli_containing_point)
	for(int starting_plg_id=0; starting_plg_id < plgs.size(); starting_plg_id++) {
		int thread_id = omp_get_thread_num();
		//cout << "Processing image " << starting_plg_id << "\n";
		for(const auto starting_polyline_id : potentially_compatible_polylines[starting_plg_id])
		{
			//cout << "Processing polyline " << starting_polyline_id << "\n";
			const PolyLineGraph2D::polyline &pl = plgs[starting_plg_id].polylines[starting_polyline_id];
			PolyLineGraph2D::polyline::pl_point plp = pl.get_start_plp();
			bool reached_end;
			int count=0;
			plp = pl.next_pl_point_by_distance(plp,pl.end,SPLIT_INTERVAL_DISTANCE,reached_end);
			while(!reached_end) {
				//cout << "PLP start: " << plp.segment_index << " - " << plp.coords << "\n";
				if(!plgmm.is_matched(starting_plg_id,starting_polyline_id,plp,pli_containing_point)) {
					count++;
					PolyLineGraph2D::plg_point starting_plgp(starting_polyline_id,plp);

					// Detect new 3D points using starting_plgp and potentially compatible polylines
					vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> cur_res = find_new_3d_points_from_compatible_polylines_starting_plgp_expandallviews(sfmd, plgs, all_fundamental_matrices, potentially_compatible_polylines, starting_plg_id, starting_plgp,plgmm,plmaps);

					for(auto &new_p3d : cur_res)
						res_t[thread_id].push_back(new_p3d);

					plp = pl.next_pl_point_by_distance(plp,pl.end,SPLIT_INTERVAL_DISTANCE,reached_end);
				} else {
					// jump to interval end
					plp = pl.next_pl_point_by_distance(pli_containing_point.end,pl.end,SPLIT_INTERVAL_DISTANCE,reached_end);
				}
			}
		}
	}

	for(int i=0; i < num_threads; i++)
		for(auto &new_p3d : res_t[i])
			res.push_back(new_p3d);

	delete[] res_t;

	return res;
}

/**
 * Given
 * vector<vector<ulong>> : for each image, vector of polyline ids
 * representing sets of polylines on all images considered potentially compatible among them
 *
 * returns new 3D points vector
**/
vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> find_new_3d_points_from_compatible_polylines_expandallviews(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const Mat** all_fundamental_matrices, const vector<set<ulong>> &potentially_compatible_polylines, PLGMatchesManager &plgmm, const vector<PolyLine2DMapSearch> &plmaps)
{
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> res;
	PolyLineGraph2D::polyline::pl_interval pli_containing_point;

	for(int starting_plg_id=0; starting_plg_id < plgs.size(); starting_plg_id++) {
		//cout << "Processing image " << starting_plg_id << "\n";
		for(const auto starting_polyline_id : potentially_compatible_polylines[starting_plg_id])
		{
			//cout << "Processing polyline " << starting_polyline_id << "\n";
			const PolyLineGraph2D::polyline &pl = plgs[starting_plg_id].polylines[starting_polyline_id];
			PolyLineGraph2D::polyline::pl_point plp = pl.get_start_plp();
			bool reached_end;
			int count=0;
			plp = pl.next_pl_point_by_distance(plp,pl.end,SPLIT_INTERVAL_DISTANCE,reached_end);
			while(!reached_end) {
				if(!plgmm.is_matched(starting_plg_id,starting_polyline_id,plp,pli_containing_point)) {
					count++;
					PolyLineGraph2D::plg_point starting_plgp(starting_polyline_id,plp);

					// Detect new 3D points using starting_plgp and potentially compatible polylines
					vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> cur_res = find_new_3d_points_from_compatible_polylines_starting_plgp_expandallviews(sfmd, plgs, all_fundamental_matrices, potentially_compatible_polylines, starting_plg_id, starting_plgp,plgmm,plmaps);

					plgmm.add_matched_3dpolyline(cur_res);

					for(auto &new_p3d : cur_res)
						res.push_back(new_p3d);

					plp = pl.next_pl_point_by_distance(plp,pl.end,SPLIT_INTERVAL_DISTANCE,reached_end);
				} else {
					// jump to interval end
					plp = pl.next_pl_point_by_distance(pli_containing_point.end,pl.end,SPLIT_INTERVAL_DISTANCE,reached_end);
				}
			}
		}
	}

	return res;
}
