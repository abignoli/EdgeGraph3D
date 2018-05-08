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


#include "plg_matching.hpp"

#include <opencv2/core/mat.hpp>

#include "SfMData.h"
#include "glm.hpp"


/**
 * - Find next point on first polyline (by distance)
 * 		if end is reached return false
 * - Find epipolar line from first PLG to the second
 * - Find next point on second polyline (by epipolar line intersection)
 * 		if end/quasiparallel is reached return false
 * - Find next point on third polyline (by epipolar line intersection)
 * 		if end/quasiparallel is reached return false
 * - Triangulate a 3D point using the computed 2D correspondences
 * 		return false if triangulation fails
 * - return true
 */
bool compatible(const SfMData &sfmd, const vector<int> &selected_2d_reprojections_ids, const tuple<PolyLineGraph2DHMapImpl,PolyLineGraph2DHMapImpl,PolyLineGraph2DHMapImpl> &plgs, const tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point> &current_plgps, const tuple<ulong,ulong,ulong> &directions, const Mat &fab, const Mat &fac, tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point> &next_plgps, std::tuple<glm::vec3, vector<glm::vec2>, vector<int>> &new_point_data, bool &fail_due_to_parallel_epipolar) {
	fail_due_to_parallel_epipolar = false;

	const PolyLineGraph2DHMapImpl &plg_a = get<0>(plgs);
	const PolyLineGraph2DHMapImpl &plg_b = get<1>(plgs);
	const PolyLineGraph2DHMapImpl &plg_c = get<2>(plgs);

	const PolyLineGraph2D::plg_point &plgp_a = get<0>(current_plgps);
	const PolyLineGraph2D::plg_point &plgp_b = get<1>(current_plgps);
	const PolyLineGraph2D::plg_point &plgp_c = get<2>(current_plgps);

	const PolyLineGraph2D::polyline &pl_a = plg_a.polylines[plgp_a.polyline_id];
	const PolyLineGraph2D::polyline &pl_b = plg_b.polylines[plgp_b.polyline_id];
	const PolyLineGraph2D::polyline &pl_c = plg_c.polylines[plgp_c.polyline_id];

	const ulong dir_a = get<0>(directions);
	const ulong dir_b = get<1>(directions);
	const ulong dir_c = get<2>(directions);

	bool reached_polyline_extreme;
	const PolyLineGraph2D::polyline::pl_point next_plp_a = pl_a.next_pl_point_by_distance(PolyLineGraph2D::polyline::pl_point(plgp_a.plp.segment_index,plgp_a.plp.coords),dir_a,PLG_FOLLOW_FIRST_IMAGE_DISTANCE,reached_polyline_extreme);
	if(reached_polyline_extreme) {
		return false;
	}

	glm::vec3 epipolar;
	bool found;
	PolyLineGraph2D::polyline::pl_point next_before_quasiparallel;

	if(!computeCorrespondEpilineSinglePoint(next_plp_a.coords, fab, epipolar, 1))
		return false;
	PolyLineGraph2D::polyline::pl_point next_plp_b;
	pl_b.next_pl_point_by_line_intersection(
			PolyLineGraph2D::polyline::pl_point(plgp_b.plp.segment_index,plgp_b.plp.coords),
			dir_b,
			epipolar,
			next_plp_b,
			fail_due_to_parallel_epipolar,
			next_before_quasiparallel,
			reached_polyline_extreme,
			found);
	if(!found) {
		return false;
	}

	if(!computeCorrespondEpilineSinglePoint(next_plp_a.coords, fac, epipolar, 1))
		return false;
	PolyLineGraph2D::polyline::pl_point next_plp_c;
	pl_c.next_pl_point_by_line_intersection(
			PolyLineGraph2D::polyline::pl_point(plgp_c.plp.segment_index,plgp_c.plp.coords),
			dir_c,
			epipolar,
			next_plp_c,
			fail_due_to_parallel_epipolar,
			next_before_quasiparallel,
			reached_polyline_extreme,
			found);
	if(!found) {
		return false;
	}

	bool valid;
	vector<glm::vec2> selected_2d_reprojections_coords;
	selected_2d_reprojections_coords.push_back(next_plp_a.coords);
	selected_2d_reprojections_coords.push_back(next_plp_b.coords);
	selected_2d_reprojections_coords.push_back(next_plp_c.coords);
	compute_3d_point(sfmd,
			selected_2d_reprojections_coords,
			selected_2d_reprojections_ids,
			new_point_data,
			valid);
	if(!valid) {
		return false;
	}

	next_plgps = make_tuple(
			PolyLineGraph2D::plg_point(plgp_a.polyline_id,next_plp_a.segment_index,next_plp_a.coords),
			PolyLineGraph2D::plg_point(plgp_b.polyline_id,next_plp_b.segment_index,next_plp_b.coords),
			PolyLineGraph2D::plg_point(plgp_c.polyline_id,next_plp_c.segment_index,next_plp_c.coords));

	return true;
}

bool compatible(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<int> &selected_2d_reprojections_ids, const vector<PolyLineGraph2DHMapImpl> &plgs, const vector<PolyLineGraph2D::plg_point> &current_plgps, const tuple<ulong,ulong,ulong> &directions, tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point> &next_plgps, std::tuple<glm::vec3, vector<glm::vec2>, vector<int>> &new_point_data, bool &fail_due_to_parallel_epipolar) {
	const tuple<PolyLineGraph2DHMapImpl,PolyLineGraph2DHMapImpl,PolyLineGraph2DHMapImpl> selected_plgs = make_tuple(plgs[selected_2d_reprojections_ids[0]],plgs[selected_2d_reprojections_ids[1]],plgs[selected_2d_reprojections_ids[2]]);
	const tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point> selected_current_plgps = make_tuple(current_plgps[selected_2d_reprojections_ids[0]],current_plgps[selected_2d_reprojections_ids[1]],current_plgps[selected_2d_reprojections_ids[2]]);;
	const Mat &fab = all_fundamental_matrices[selected_2d_reprojections_ids[0]][selected_2d_reprojections_ids[1]];
	const Mat &fac = all_fundamental_matrices[selected_2d_reprojections_ids[0]][selected_2d_reprojections_ids[2]];
	return compatible(sfmd, selected_2d_reprojections_ids, selected_plgs, selected_current_plgps, directions, fab, fac, next_plgps, new_point_data, fail_due_to_parallel_epipolar);
}

bool find_direction_given_first_extreme(const SfMData &sfmd, const vector<int> &selected_2d_reprojections_ids, const tuple<PolyLineGraph2DHMapImpl,PolyLineGraph2DHMapImpl,PolyLineGraph2DHMapImpl> &plgs, const tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point> &current_plgps, const Mat &fab, const Mat &fac, const ulong first_direction,tuple<ulong,ulong,ulong> &valid_direction, vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> &valid_points, bool &fail_due_to_parallel_epipolar)
{
	const PolyLineGraph2DHMapImpl &plg_a = get<0>(plgs);
	const PolyLineGraph2DHMapImpl &plg_b = get<1>(plgs);
	const PolyLineGraph2DHMapImpl &plg_c = get<2>(plgs);

	const PolyLineGraph2D::plg_point &plgp_a = get<0>(current_plgps);
	const PolyLineGraph2D::plg_point &plgp_b = get<1>(current_plgps);
	const PolyLineGraph2D::plg_point &plgp_c = get<2>(current_plgps);

	const PolyLineGraph2D::polyline &pl_a = plg_a.polylines[plgp_a.polyline_id];
	const PolyLineGraph2D::polyline &pl_b = plg_b.polylines[plgp_b.polyline_id];
	const PolyLineGraph2D::polyline &pl_c = plg_c.polylines[plgp_c.polyline_id];

	vector<tuple<ulong,ulong,ulong>> potential_directions;

	potential_directions.push_back(make_tuple(first_direction,pl_b.start,pl_c.start));
	potential_directions.push_back(make_tuple(first_direction,pl_b.start,pl_c.end));
	potential_directions.push_back(make_tuple(first_direction,pl_b.end,pl_c.start));
	potential_directions.push_back(make_tuple(first_direction,pl_b.end,pl_c.end));

	vector<bool> valid_directions;
	vector<vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>>> triangulated_points;
	vector<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>> last_plgps;
	for(const auto &pd : potential_directions) {
		valid_directions.push_back(true);
		last_plgps.push_back(tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>(current_plgps));
		triangulated_points.push_back(vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>>());
	}
	int amount_of_valid = valid_directions.size();

	bool first_it = true;
	while(amount_of_valid > 1) {
		for(int i=0; i < potential_directions.size(); i++) {
			if(valid_directions[i]) {
				const tuple<ulong,ulong,ulong> &pd = potential_directions[i];
				std::tuple<glm::vec3, vector<glm::vec2>, vector<int>> new_point_data;
				tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point> next_plgps;
				if(compatible(sfmd,selected_2d_reprojections_ids,plgs,last_plgps[i],pd,fab,fac,next_plgps,new_point_data, fail_due_to_parallel_epipolar)) {
					last_plgps[i] = next_plgps;
					triangulated_points[i].push_back(make_pair(next_plgps,new_point_data));
				} else {
					valid_directions[i] = false;
					amount_of_valid--;
				}
			}
		}
		first_it = false;
	}

	if(amount_of_valid == 0)
		// Fail
		return false;
	else {
		for(int i=0; i < potential_directions.size(); i++)
			if(valid_directions[i]) {
				valid_direction = potential_directions[i];
				valid_points = triangulated_points[i];
			}
		return true;
	}
}

void find_directions(const SfMData &sfmd,  const vector<int> &selected_2d_reprojections_ids, const tuple<PolyLineGraph2DHMapImpl,PolyLineGraph2DHMapImpl,PolyLineGraph2DHMapImpl> &plgs, const tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point> &current_plgps, const Mat &fab, const Mat &fac, tuple<ulong,ulong,ulong> &direction1, bool &direction1_valid, vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> &valid_points_direction1, tuple<ulong,ulong,ulong> &direction2, bool &direction2_valid, vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> &valid_points_direction2, bool &fail_due_to_parallel_epipolar) {
	direction1_valid = false;
	direction2_valid = false;

	const PolyLineGraph2DHMapImpl &plg_a = get<0>(plgs);
	const PolyLineGraph2DHMapImpl &plg_b = get<1>(plgs);
	const PolyLineGraph2DHMapImpl &plg_c = get<2>(plgs);

	const PolyLineGraph2D::plg_point &plgp_a = get<0>(current_plgps);
	const PolyLineGraph2D::plg_point &plgp_b = get<1>(current_plgps);
	const PolyLineGraph2D::plg_point &plgp_c = get<2>(current_plgps);

	const PolyLineGraph2D::polyline &pl_a = plg_a.polylines[plgp_a.polyline_id];
	const PolyLineGraph2D::polyline &pl_b = plg_b.polylines[plgp_b.polyline_id];
	const PolyLineGraph2D::polyline &pl_c = plg_c.polylines[plgp_c.polyline_id];

	vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> valid_points_towards_first_start;

	bool valid_towards_start = find_direction_given_first_extreme(sfmd,selected_2d_reprojections_ids,plgs,current_plgps,fab,fac,pl_a.start,direction1, valid_points_towards_first_start, fail_due_to_parallel_epipolar);
	if(valid_towards_start) {
		// It was possible to find a valid directions configuration towards the start of the polyline on the first image
		direction1_valid = true;
		valid_points_direction1 = valid_points_towards_first_start;

		// Get opposite direction
		direction2 = make_tuple(
				pl_a.start == get<0>(direction1) ? pl_a.end : pl_a.start,
				pl_b.start == get<1>(direction1) ? pl_b.end : pl_b.start,
				pl_c.start == get<2>(direction1) ? pl_c.end : pl_c.start);

		// Test opposite direction
		std::tuple<glm::vec3, vector<glm::vec2>, vector<int>> new_point_data;
		tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point> next_plgps;
		if(compatible(sfmd,selected_2d_reprojections_ids,plgs,current_plgps,direction2,fab,fac,next_plgps,new_point_data, fail_due_to_parallel_epipolar)) {
			// Other direction is also compatible
			direction2_valid = true;
			vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> valid_points_opposite;
			valid_points_opposite.push_back(make_pair(next_plgps,new_point_data));
			valid_points_direction2 = valid_points_opposite;
		}

	} else {
		// It was not possible to find a valid directions configuration towards the start of the polyline on the first image

		vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> valid_points_towards_first_end;

		bool valid_towards_end = find_direction_given_first_extreme(sfmd,selected_2d_reprojections_ids,plgs,current_plgps,fab,fac,pl_a.end,direction1, valid_points_towards_first_end, fail_due_to_parallel_epipolar);

		if(valid_towards_end) {
			// A valid direction has been found towards the first end
			direction1_valid = true;
			valid_points_direction1 = valid_points_towards_first_end;

			// Get opposite direction
			direction2 = make_tuple(
					pl_a.start == get<0>(direction1) ? pl_a.end : pl_a.start,
					pl_b.start == get<1>(direction1) ? pl_b.end : pl_b.start,
					pl_c.start == get<2>(direction1) ? pl_c.end : pl_c.start);
		}
	}
}

void find_directions(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs, const std::tuple<glm::vec3,vector<PolyLineGraph2D::plg_point>,vector<int>> &matches, const int selected_indexes[], tuple<ulong,ulong,ulong> &direction1, bool &direction1_valid, vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> &valid_points_direction1, tuple<ulong,ulong,ulong> &direction2, bool &direction2_valid, vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> &valid_points_direction2, bool &fail_due_to_parallel_epipolar) {
	const vector<PolyLineGraph2D::plg_point> &matches_data = get<1>(matches);
	const vector<int> &matched_ids = get<2>(matches);
	vector<int> selected_2d_reprojections_ids;
	selected_2d_reprojections_ids.push_back(matched_ids[selected_indexes[0]]);
	selected_2d_reprojections_ids.push_back(matched_ids[selected_indexes[1]]);
	selected_2d_reprojections_ids.push_back(matched_ids[selected_indexes[2]]);
	const tuple<PolyLineGraph2DHMapImpl,PolyLineGraph2DHMapImpl,PolyLineGraph2DHMapImpl> selected_plgs = make_tuple(plgs[selected_2d_reprojections_ids[0]],plgs[selected_2d_reprojections_ids[1]],plgs[selected_2d_reprojections_ids[2]]);
	const tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point> selected_current_plgps = make_tuple(matches_data[selected_indexes[0]],matches_data[selected_indexes[1]],matches_data[selected_indexes[2]]);;
	const Mat &fab = all_fundamental_matrices[selected_2d_reprojections_ids[0]][selected_2d_reprojections_ids[1]];
	const Mat &fac = all_fundamental_matrices[selected_2d_reprojections_ids[0]][selected_2d_reprojections_ids[2]];
	find_directions(sfmd, selected_2d_reprojections_ids, selected_plgs, selected_current_plgps, fab, fac, direction1, direction1_valid, valid_points_direction1, direction2, direction2_valid, valid_points_direction2, fail_due_to_parallel_epipolar);
}

/**
 * Find directions only using first and last view
 */
void find_directions_3view_firstlast(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs, const std::tuple<glm::vec3,vector<PolyLineGraph2D::plg_point>,vector<int>> &matches, tuple<ulong,ulong,ulong> &direction1, bool &direction1_valid, vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> &valid_points_direction1, tuple<ulong,ulong,ulong> &direction2, bool &direction2_valid, vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> &valid_points_direction2, bool &fail_due_to_parallel_epipolar) {
	const vector<PolyLineGraph2D::plg_point> &matches_data = get<1>(matches);
	const vector<int> &matched_ids = get<2>(matches);
	int selected_indexes[3];
	const int amount_of_matches = matched_ids.size();

	selected_indexes[0] = 0;
	selected_indexes[1] = amount_of_matches/2;
	selected_indexes[2] = amount_of_matches-1;

	vector<int> selected_2d_reprojections_ids;
	selected_2d_reprojections_ids.push_back(matched_ids[selected_indexes[0]]);
	selected_2d_reprojections_ids.push_back(matched_ids[selected_indexes[1]]);
	selected_2d_reprojections_ids.push_back(matched_ids[selected_indexes[2]]);

	find_directions(sfmd,all_fundamental_matrices, plgs,matches, selected_indexes, direction1, direction1_valid, valid_points_direction1, direction2, direction2_valid, valid_points_direction2, fail_due_to_parallel_epipolar);
}

/**
 * Find directions only using first and last view
 */
void find_directions_3view_firstlast(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs, const std::tuple<glm::vec3,vector<PolyLineGraph2D::plg_point>,vector<int>> &matches, tuple<ulong,ulong,ulong> &direction1, bool &direction1_valid, vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> &valid_points_direction1, tuple<ulong,ulong,ulong> &direction2, bool &direction2_valid, vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> &valid_points_direction2, bool &fail_due_to_parallel_epipolar, vector<int> &selected_2d_reprojections_ids) {
	const vector<PolyLineGraph2D::plg_point> &matches_data = get<1>(matches);
	const vector<int> &matched_ids = get<2>(matches);
	int selected_indexes[3];
	const int amount_of_matches = matched_ids.size();

	selected_indexes[0] = 0;
	selected_indexes[1] = amount_of_matches/2;
	selected_indexes[2] = amount_of_matches-1;

	selected_2d_reprojections_ids.push_back(matched_ids[selected_indexes[0]]);
	selected_2d_reprojections_ids.push_back(matched_ids[selected_indexes[1]]);
	selected_2d_reprojections_ids.push_back(matched_ids[selected_indexes[2]]);

	find_directions(sfmd,all_fundamental_matrices, plgs,matches, selected_indexes, direction1, direction1_valid, valid_points_direction1, direction2, direction2_valid, valid_points_direction2, fail_due_to_parallel_epipolar);
}

/**
 * Find directions only using first and last view
 */
void find_directions_3view_firstlast(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs, const std::tuple<glm::vec3,vector<PolyLineGraph2D::plg_point>,vector<int>> &matches, vector<ulong> &direction1, bool &direction1_valid, vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &valid_points_direction1, vector<ulong> &direction2, bool &direction2_valid, vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &valid_points_direction2, bool &fail_due_to_parallel_epipolar) {
	tuple<ulong,ulong,ulong> direction1_t,direction2_t;
	vector<int> selected_2d_reprojections_ids;

	vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> valid_points_direction1_t, valid_points_direction2_t;
	find_directions_3view_firstlast(sfmd, all_fundamental_matrices, plgs, matches, direction1_t, direction1_valid, valid_points_direction1_t, direction2_t, direction2_valid, valid_points_direction2_t, fail_due_to_parallel_epipolar,selected_2d_reprojections_ids);
	if(direction1_valid) {
		direction1 = vector<ulong>(plgs.size());

		direction1[selected_2d_reprojections_ids[0]] = get<0>(direction1_t);
		direction1[selected_2d_reprojections_ids[1]] = get<1>(direction1_t);
		direction1[selected_2d_reprojections_ids[2]] = get<2>(direction1_t);

		valid_points_direction1.clear();
		for(auto &vpt : valid_points_direction1_t) {
			vector<PolyLineGraph2D::plg_point> selected_plgps;
			selected_plgps.push_back(get<0>(vpt.first));
			selected_plgps.push_back(get<1>(vpt.first));
			selected_plgps.push_back(get<2>(vpt.first));

			valid_points_direction1.push_back(make_tuple(
					get<0>(vpt.second),
					selected_plgps,
					get<2>(vpt.second)));
		}

		direction2 = vector<ulong>(plgs.size());
		direction2[selected_2d_reprojections_ids[0]] = get<0>(direction2_t);
		direction2[selected_2d_reprojections_ids[1]] = get<1>(direction2_t);
		direction2[selected_2d_reprojections_ids[2]] = get<2>(direction2_t);
		if(direction2_valid) {
			valid_points_direction2.clear();
			for(auto &vpt : valid_points_direction2_t) {
				vector<PolyLineGraph2D::plg_point> selected_plgps;
				selected_plgps.push_back(get<0>(vpt.first));
				selected_plgps.push_back(get<1>(vpt.first));
				selected_plgps.push_back(get<2>(vpt.first));

				valid_points_direction2.push_back(make_tuple(
						get<0>(vpt.second),
						selected_plgps,
						get<2>(vpt.second)));
			}
		}
	}
}

/**
 * Find directions only using first and last view. If that fails, try matching from last to first
 */
void find_directions_3view_firstlast_lastfirst(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs, const std::tuple<glm::vec3,vector<PolyLineGraph2D::plg_point>,vector<int>> &matches, tuple<ulong,ulong,ulong> &direction1, bool &direction1_valid, vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> &valid_points_direction1, tuple<ulong,ulong,ulong> &direction2, bool &direction2_valid, vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> &valid_points_direction2, bool &fail_due_to_parallel_epipolar) {
	const vector<PolyLineGraph2D::plg_point> &matches_data = get<1>(matches);
	const vector<int> &matched_ids = get<2>(matches);
	const int amount_of_matches = matched_ids.size();

	int selected_indexes[3];

	selected_indexes[0] = 0;
	selected_indexes[1] = amount_of_matches/2;
	selected_indexes[2] = amount_of_matches-1;

	find_directions(sfmd,all_fundamental_matrices, plgs,matches, selected_indexes, direction1, direction1_valid, valid_points_direction1, direction2, direction2_valid, valid_points_direction2, fail_due_to_parallel_epipolar);

	if(!direction2_valid) {
		// Failed to follow in both direction
		if(fail_due_to_parallel_epipolar) {
			int tmp = selected_indexes[0];
			selected_indexes[0] = selected_indexes[2];
			selected_indexes[2] = tmp;

			if(!direction1_valid) {
				// Could not detect directions due to epipolar paralleliness in either direction
				find_directions(sfmd,all_fundamental_matrices, plgs,matches, selected_indexes, direction1, direction1_valid, valid_points_direction1, direction2, direction2_valid, valid_points_direction2, fail_due_to_parallel_epipolar);
			} else {
				// Could not detect directions due to epipolar paralleliness just in second direction

				vector<int> selected_2d_reprojections_ids;
				selected_2d_reprojections_ids.push_back(matched_ids[selected_indexes[0]]);
				selected_2d_reprojections_ids.push_back(matched_ids[selected_indexes[1]]);
				selected_2d_reprojections_ids.push_back(matched_ids[selected_indexes[2]]);
				const tuple<PolyLineGraph2DHMapImpl,PolyLineGraph2DHMapImpl,PolyLineGraph2DHMapImpl> selected_plgs = make_tuple(plgs[selected_2d_reprojections_ids[0]],plgs[selected_2d_reprojections_ids[1]],plgs[selected_2d_reprojections_ids[2]]);
				const tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point> selected_current_plgps = make_tuple(matches_data[selected_indexes[0]],matches_data[selected_indexes[1]],matches_data[selected_indexes[2]]);;
				const Mat &fab = all_fundamental_matrices[selected_2d_reprojections_ids[0]][selected_2d_reprojections_ids[1]];
				const Mat &fac = all_fundamental_matrices[selected_2d_reprojections_ids[0]][selected_2d_reprojections_ids[2]];
				const PolyLineGraph2DHMapImpl &plg_a = get<0>(selected_plgs);
				const PolyLineGraph2DHMapImpl &plg_b = get<1>(selected_plgs);
				const PolyLineGraph2DHMapImpl &plg_c = get<2>(selected_plgs);

				const PolyLineGraph2D::plg_point &plgp_a = get<0>(selected_current_plgps);
				const PolyLineGraph2D::plg_point &plgp_b = get<1>(selected_current_plgps);
				const PolyLineGraph2D::plg_point &plgp_c = get<2>(selected_current_plgps);

				const PolyLineGraph2D::polyline &pl_a = plg_a.polylines[plgp_a.polyline_id];
				const PolyLineGraph2D::polyline &pl_b = plg_b.polylines[plgp_b.polyline_id];
				const PolyLineGraph2D::polyline &pl_c = plg_c.polylines[plgp_c.polyline_id];

				// Test the opposite to direction 1, but on last-first matching

				// Get opposite direction (reversed, now first is last i.e. what now is plg_a was plg_c in the earlier computation)
				tuple<ulong,ulong,ulong> opposite_direction = make_tuple(
						pl_a.start == get<2>(direction1) ? pl_a.end : pl_a.start,
						pl_b.start == get<1>(direction1) ? pl_b.end : pl_b.start,
						pl_c.start == get<0>(direction1) ? pl_c.end : pl_c.start);

				// Test opposite direction
				std::tuple<glm::vec3, vector<glm::vec2>, vector<int>> new_point_data;
				tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point> next_plgps;
				if(compatible(sfmd,selected_2d_reprojections_ids,selected_plgs,selected_current_plgps,opposite_direction,fab,fac,next_plgps,new_point_data, fail_due_to_parallel_epipolar)) {
					// cout << "----> Second direction identified thanks to last-first matching\n";

					//tmp_lastfirst_newdir_counter++;

					// Other direction is also compatible
					direction2 = make_tuple(get<2>(opposite_direction),get<1>(opposite_direction),get<0>(opposite_direction));
					direction2_valid = true;
					vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> valid_points_opposite;
					tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point> reverse_next_plgps = make_tuple(get<0>(next_plgps),get<1>(next_plgps),get<2>(next_plgps));
					std::reverse(get<1>(new_point_data).begin(),get<1>(new_point_data).end());
					std::reverse(get<2>(new_point_data).begin(),get<2>(new_point_data).end());
					valid_points_opposite.push_back(make_pair(reverse_next_plgps,new_point_data));
					valid_points_direction2 = valid_points_opposite;
				}
			}
		}
	}
}

#if defined(FOLLOW_DIRECTION_COMPATIBLE_FORCE_EXTREME_MATCHING)
/**
 * - Find next point on first polyline (by distance)
 * 		if end is reached return false
 * - Find epipolar line from first PLG to the second
 * - Find next point on second polyline (by epipolar line intersection)
 * 		if end/quasiparallel is reached return false
 * - Find next point on third polyline (by epipolar line intersection)
 * 		if end/quasiparallel is reached return false
 * - Triangulate a 3D point using the computed 2D correspondences
 * 		return false if triangulation fails
 * - return true
*/
bool compatible(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const vector<ulong> &directions, const Mat** all_fundamental_matrices, const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &current_plgp, std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &new_point_data, bool &fail_due_to_parallel_epipolar) {
	//const int starting_plg_index=0;
	for(int starting_plg_index=0; starting_plg_index < get<2>(current_plgp).size(); starting_plg_index++) {
		// cout << "Checking next point compatibility\n";

		const vector<PolyLineGraph2D::plg_point> &current_plgps = get<1>(current_plgp);
		const vector<int> &current_plgps_ids = get<2>(current_plgp);

	//	unsigned int amount_of_potential_correspondences = 0;

		const int starting_plg_id = current_plgps_ids[starting_plg_index];
		vector<int> selected_2d_reprojections_ids;
		vector<glm::vec2> selected_2d_reprojections_coords;
		vector<PolyLineGraph2D::plg_point> selected_plgps;

		// Get next plgp on starting plg (by distance)
		bool reached_polyline_extreme;
		const PolyLineGraph2DHMapImpl &plg_starting = plgs[starting_plg_id];
		const PolyLineGraph2D::plg_point &plgp_starting = current_plgps[starting_plg_index];
		const PolyLineGraph2D::polyline &pl_starting = plg_starting.polylines[plgp_starting.polyline_id];
		const PolyLineGraph2D::polyline::pl_point next_starting_plg_plp = pl_starting.next_pl_point_by_distance(PolyLineGraph2D::polyline::pl_point(plgp_starting.segment_index,plgp_starting.coords),directions[starting_plg_id],PLG_FOLLOW_FIRST_IMAGE_DISTANCE,reached_polyline_extreme);
		// cout << "Searching first 2D Correspondence found by distance on PLG " << starting_plg_id << " from " << plgp_starting.coords << "\n";
		if(reached_polyline_extreme) {
			// cout << "First PLG (" << starting_plg_id << ") reached extreme. Trying to match extreme\n";
			//continue;
		}
		// cout << "First 2D Correspondence found by distance on PLG " << starting_plg_id << " at " << next_starting_plg_plp.coords << "\n";


	//	amount_of_potential_correspondences++;
		selected_plgps.push_back(PolyLineGraph2D::plg_point(plgp_starting.polyline_id,next_starting_plg_plp.segment_index,next_starting_plg_plp.coords));
		selected_2d_reprojections_ids.push_back(starting_plg_id);
		selected_2d_reprojections_coords.push_back(next_starting_plg_plp.coords);

		glm::vec3 epipolar;
		bool found;
		bool bounded_distance_violated;
		bool need_restart = false;
		PolyLineGraph2D::polyline::pl_point next_before_quasiparallel;

		// Get next plgp on other plgs (by line intersection)
		for(int i=0; i < current_plgps_ids.size(); i++)
			if(i != starting_plg_index) {

				const int cur_plg_id = current_plgps_ids[i];
				const PolyLineGraph2D::plg_point &cur_plgp = current_plgps[i];
				const PolyLineGraph2D::polyline &cur_pl = plgs[cur_plg_id].polylines[cur_plgp.polyline_id];

				// cout << "Searching 2D Correspondence found by epipolar-intersection on PLG " << cur_plg_id << " from " << cur_plgp.coords << "\n";

				computeCorrespondEpilineSinglePoint(next_starting_plg_plp.coords, all_fundamental_matrices[starting_plg_id][cur_plg_id], epipolar, 1);

				PolyLineGraph2D::polyline::pl_point next_plp;
/*				cur_pl.next_pl_point_by_line_intersection(
						PolyLineGraph2D::polyline::pl_point(cur_plgp.segment_index,cur_plgp.coords),
						directions[cur_plg_id],
						epipolar,
						next_plp,
						fail_due_to_parallel_epipolar,
						next_before_quasiparallel,
						reached_polyline_extreme,
						found);*/
				cur_pl.next_pl_point_by_line_intersection_bounded_distance(
						PolyLineGraph2D::polyline::pl_point(cur_plgp.segment_index,cur_plgp.coords),
						directions[cur_plg_id],
						epipolar,
						PLG_FOLLOW_CORRESPONDENCE_IMAGE_DISTANCE_MIN,
						PLG_FOLLOW_CORRESPONDENCE_IMAGE_DISTANCE_MAX,
						next_plp,
						fail_due_to_parallel_epipolar,
						next_before_quasiparallel,
						reached_polyline_extreme,
						bounded_distance_violated,
						found);
				if(found) {
	//				amount_of_potential_correspondences++;
					selected_plgps.push_back(PolyLineGraph2D::plg_point(cur_plgp.polyline_id,next_plp.segment_index,next_plp.coords));
					selected_2d_reprojections_ids.push_back(cur_plg_id);
					selected_2d_reprojections_coords.push_back(next_plp.coords);
					// cout << "2D Correspondence found by epipolar-intersection on PLG " << cur_plg_id << " at " << next_plp.coords << "\n";
				} else {
					// cout << "Could not find 2D Correspondence on PLG " << cur_plg_id << ": ";
					if(reached_polyline_extreme) {
						// cout << "reached polyline extreme. ";
						if(i > starting_plg_index) {
							starting_plg_index = i-1;
							need_restart = true;
							// cout << "Re-starting from this PLG\n";
						} else
							// cout << "Not re-starting.\n";
						break;
					} else if (fail_due_to_parallel_epipolar) {
						// cout << "found quasi-parallel segment before or while reaching en epipolar intersection\n";
					}
				}
			}

		if(need_restart)
			continue;

		bool valid;

		glm::vec3 new_3d_coords;

		if(selected_2d_reprojections_ids.size() < PLG_MATCHING_TRIANGULATION_MINIMUM_AMOUNT_OF_POINTS) {
			// cout << "Not enough 2D Correspondences to continue following PLG\n";
			continue;
			//return false;
		}

		compute_3d_point_coords(sfmd,
				selected_2d_reprojections_coords,
				selected_2d_reprojections_ids,
				new_3d_coords,
				valid);
		if(!valid) {
			vector<bool> selected;
			vector<glm::vec2> actually_selected_2d_reprojections_coords;
			vector<int> actually_selected_2d_reprojections_ids;
			compute_3d_point_coords_combinations(sfmd,
				selected_2d_reprojections_coords,
				selected_2d_reprojections_ids,
				PLG_MATCHING_TRIANGULATION_MINIMUM_AMOUNT_OF_POINTS,
				actually_selected_2d_reprojections_coords,
				actually_selected_2d_reprojections_ids,
				selected,
				new_3d_coords,
				valid);

			if(valid) {
				selected_2d_reprojections_coords = actually_selected_2d_reprojections_coords;
				selected_2d_reprojections_ids = actually_selected_2d_reprojections_ids;

				vector<PolyLineGraph2D::plg_point> actually_selected_plgps;
				for(int i=0; i < selected_plgps.size(); i++)
					if(selected[i])
						actually_selected_plgps.push_back(selected_plgps[i]);
				selected_plgps = actually_selected_plgps;
			}
		}

		if(valid) {
			get<0>(new_point_data) = new_3d_coords;
			get<1>(new_point_data) = selected_plgps;
			get<2>(new_point_data) = selected_2d_reprojections_ids;

			return true;
		}
	}

	// cout << "Could not triangulate 2D correspondences given by directions\n";
	return false;
}
#else
/*
*
 * - Find next point on first polyline (by distance)
 * 		if end is reached return false
 * - Find epipolar line from first PLG to the second
 * - Find next point on second polyline (by epipolar line intersection)
 * 		if end/quasiparallel is reached return false
 * - Find next point on third polyline (by epipolar line intersection)
 * 		if end/quasiparallel is reached return false
 * - Triangulate a 3D point using the computed 2D correspondences
 * 		return false if triangulation fails
 * - return true
*/

bool compatible(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const vector<ulong> &directions, const Mat** all_fundamental_matrices, const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &current_plgp, std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &new_point_data, bool &fail_due_to_parallel_epipolar) {
	//const int starting_plg_index=0;
	for(int starting_plg_index=0; starting_plg_index < get<2>(current_plgp).size(); starting_plg_index++) {
		// cout << "Checking next point compatibility\n";

		const vector<PolyLineGraph2D::plg_point> &current_plgps = get<1>(current_plgp);
		const vector<int> &current_plgps_ids = get<2>(current_plgp);

	//	unsigned int amount_of_potential_correspondences = 0;

		const int starting_plg_id = current_plgps_ids[starting_plg_index];
		vector<int> selected_2d_reprojections_ids;
		vector<glm::vec2> selected_2d_reprojections_coords;
		vector<PolyLineGraph2D::plg_point> selected_plgps;

		// Get next plgp on starting plg (by distance)
		bool reached_polyline_extreme;
		const PolyLineGraph2DHMapImpl &plg_starting = plgs[starting_plg_id];
		const PolyLineGraph2D::plg_point &plgp_starting = current_plgps[starting_plg_index];
		const PolyLineGraph2D::polyline &pl_starting = plg_starting.polylines[plgp_starting.polyline_id];
		const PolyLineGraph2D::polyline::pl_point next_starting_plg_plp = pl_starting.next_pl_point_by_distance(PolyLineGraph2D::polyline::pl_point(plgp_starting.plp.segment_index,plgp_starting.plp.coords),directions[starting_plg_id],PLG_FOLLOW_FIRST_IMAGE_DISTANCE,reached_polyline_extreme);
		// cout << "Searching first 2D Correspondence found by distance on PLG " << starting_plg_id << " from " << plgp_starting.coords << "\n";
		if(reached_polyline_extreme) {
			// cout << "First PLG (" << starting_plg_id << ") reached extreme\n";
			continue;
			//return false;
		}
		//cout << "First 2D Correspondence found by distance on PLG " << starting_plg_id << " at " << next_starting_plg_plp.coords << "\n";


	//	amount_of_potential_correspondences++;
		selected_plgps.push_back(PolyLineGraph2D::plg_point(plgp_starting.polyline_id,next_starting_plg_plp.segment_index,next_starting_plg_plp.coords));
		selected_2d_reprojections_ids.push_back(starting_plg_id);
		selected_2d_reprojections_coords.push_back(next_starting_plg_plp.coords);

		glm::vec3 epipolar;
		bool found;
		bool bounded_distance_violated;
		PolyLineGraph2D::polyline::pl_point next_before_quasiparallel;

		// Get next plgp on other plgs (by line intersection)
		for(int i=0; i < current_plgps_ids.size(); i++)
			if(i != starting_plg_index) {

				const int cur_plg_id = current_plgps_ids[i];
				const PolyLineGraph2D::plg_point &cur_plgp = current_plgps[i];
				const PolyLineGraph2D::polyline &cur_pl = plgs[cur_plg_id].polylines[cur_plgp.polyline_id];

				//cout << "Searching 2D Correspondence found by epipolar-intersection on PLG " << cur_plg_id << " from " << cur_plgp.coords << "\n";

				if(!computeCorrespondEpilineSinglePoint(next_starting_plg_plp.coords, all_fundamental_matrices[starting_plg_id][cur_plg_id], epipolar, 1))
					continue;

				PolyLineGraph2D::polyline::pl_point next_plp;

				cur_pl.next_pl_point_by_line_intersection_bounded_distance(
						PolyLineGraph2D::polyline::pl_point(cur_plgp.plp.segment_index,cur_plgp.plp.coords),
						directions[cur_plg_id],
						epipolar,
						PLG_FOLLOW_CORRESPONDENCE_IMAGE_DISTANCE_MIN,
						PLG_FOLLOW_CORRESPONDENCE_IMAGE_DISTANCE_MAX,
						next_plp,
						fail_due_to_parallel_epipolar,
						next_before_quasiparallel,
						reached_polyline_extreme,
						bounded_distance_violated,
						found);
				if(found) {
	//				amount_of_potential_correspondences++;
					selected_plgps.push_back(PolyLineGraph2D::plg_point(cur_plgp.polyline_id,next_plp.segment_index,next_plp.coords));
					selected_2d_reprojections_ids.push_back(cur_plg_id);
					selected_2d_reprojections_coords.push_back(next_plp.coords);
				}
			}

		bool valid;

		glm::vec3 new_3d_coords;

		if(selected_2d_reprojections_ids.size() < PLG_MATCHING_TRIANGULATION_MINIMUM_AMOUNT_OF_POINTS)
			continue;

		compute_3d_point_coords(sfmd,
				selected_2d_reprojections_coords,
				selected_2d_reprojections_ids,
				new_3d_coords,
				valid);
		if(!valid) {
			vector<bool> selected;
			vector<glm::vec2> actually_selected_2d_reprojections_coords;
			vector<int> actually_selected_2d_reprojections_ids;
			compute_3d_point_coords_combinations(sfmd,
				selected_2d_reprojections_coords,
				selected_2d_reprojections_ids,
				PLG_MATCHING_TRIANGULATION_MINIMUM_AMOUNT_OF_POINTS,
				actually_selected_2d_reprojections_coords,
				actually_selected_2d_reprojections_ids,
				selected,
				new_3d_coords,
				valid);

			if(valid) {
				selected_2d_reprojections_coords = actually_selected_2d_reprojections_coords;
				actually_selected_2d_reprojections_ids.clear();

				vector<PolyLineGraph2D::plg_point> actually_selected_plgps;
				for(int i=0; i < selected_plgps.size(); i++)
					if(selected[i]) {
						actually_selected_plgps.push_back(selected_plgps[i]);
						actually_selected_2d_reprojections_ids.push_back(selected_2d_reprojections_ids[i]);
					}
				selected_plgps = actually_selected_plgps;
				selected_2d_reprojections_ids = actually_selected_2d_reprojections_ids;
			}
		}

		if(valid) {
			get<0>(new_point_data) = new_3d_coords;
			get<1>(new_point_data) = selected_plgps;
			get<2>(new_point_data) = selected_2d_reprojections_ids;

			return true;
		}
	}

	return false;
}
#endif

/**
 * Follow given directions on all plgs until valid.
 */
void follow_direction(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const vector<ulong> &directions, const Mat** all_fundamental_matrices, vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &valid_points, bool &fail_due_to_parallel_epipolar) {
	std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> new_point_data;
	while(compatible(sfmd, plgs, directions, all_fundamental_matrices, valid_points[valid_points.size()-1], new_point_data, fail_due_to_parallel_epipolar))
		valid_points.push_back(std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>(new_point_data));
}

void follow_direction_vector_start(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const vector<ulong> &directions, const Mat** all_fundamental_matrices, vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &valid_points, bool &fail_due_to_parallel_epipolar) {
	std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> new_point_data;
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> new_valid_points;

	if(compatible(sfmd, plgs, directions, all_fundamental_matrices, valid_points[0], new_point_data, fail_due_to_parallel_epipolar)) {
		new_valid_points.push_back(std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>(new_point_data));

		while(compatible(sfmd, plgs, directions, all_fundamental_matrices, new_valid_points[new_valid_points.size()-1], new_point_data, fail_due_to_parallel_epipolar))
			new_valid_points.push_back(std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>(new_point_data));

		vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> res;
		for(int i=new_valid_points.size()-1;i>=0;i--)
			res.push_back(new_valid_points[i]);
		for(int i=0;i<valid_points.size();i++)
			res.push_back(valid_points[i]);

		valid_points = res;
	}
}

void follow_direction_vector_end(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const vector<ulong> &directions, const Mat** all_fundamental_matrices, vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &valid_points, bool &fail_due_to_parallel_epipolar) {
	std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> new_point_data;
	while(compatible(sfmd, plgs, directions, all_fundamental_matrices, valid_points[valid_points.size()-1], new_point_data, fail_due_to_parallel_epipolar))
		valid_points.push_back(std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>(new_point_data));
}

void get_plgp_by_epipolar_intersection_from_known_point(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const int current_plg_id, const PolyLineGraph2D::plg_point &current_plgp, const ulong direction, const Mat** all_fundamental_matrices, const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &known_point, PolyLineGraph2D::polyline::pl_point &next_plp, bool &valid) {
	const int starting_index = 0;
	const int starting_plg_id = get<2>(known_point)[starting_index];
	const glm::vec2 starting_coords = get<1>(known_point)[starting_index].plp.coords;

	glm::vec3 epipolar;
	if(!computeCorrespondEpilineSinglePoint(starting_coords, all_fundamental_matrices[starting_plg_id][current_plg_id], epipolar, 1)) {
		valid=false;
		return;
	}

	bool fail_due_to_parallel_epipolar, reached_polyline_extreme;
	PolyLineGraph2D::polyline::pl_point next_before_quasiparallel;
	plgs[current_plg_id].polylines[current_plgp.polyline_id].next_pl_point_by_line_intersection(
			PolyLineGraph2D::polyline::pl_point(current_plgp.plp.segment_index,current_plgp.plp.coords),
			direction,
			epipolar,
			next_plp,
			fail_due_to_parallel_epipolar,
			next_before_quasiparallel,
			reached_polyline_extreme,
			valid);
}

// direction on specified plg, pl compatible with given 3D points? Return amount of compatible points. Doesn't update compatible 3D points
bool compatible_direction_noupdate(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const int current_plg_id, const PolyLineGraph2D::plg_point &current_plgp, const ulong direction, const Mat** all_fundamental_matrices, const vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &valid_points, vector<std::pair<glm::vec3, PolyLineGraph2D::plg_point>> &new_points_data_to_add, bool &all_points_compatible) {
	new_points_data_to_add = vector<std::pair<glm::vec3, PolyLineGraph2D::plg_point>>();

	if(valid_points.size() == 0)
		return 0;

	bool valid;
	vector<PolyLineGraph2D::polyline::pl_point> next_plps;
	vector<glm::vec3> triangulated_points;

	PolyLineGraph2D::plg_point actual_current_plgp = current_plgp;

	for(int i=0; i < valid_points.size(); i++) {
		const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &cur_pt = valid_points[i];

		PolyLineGraph2D::polyline::pl_point next_plp;
		get_plgp_by_epipolar_intersection_from_known_point(sfmd, plgs, current_plg_id, actual_current_plgp, direction, all_fundamental_matrices, cur_pt, next_plp, valid);

		if(!valid)
			break; // fail
		else {
			// check compatibility between new plp and known 3D point observations
			glm::vec3 triangulated_point;
			em_add_new_observation_to_3Dpositions(sfmd,cur_pt,next_plp.coords,current_plg_id,triangulated_point,valid);

			if(valid) {
				triangulated_points.push_back(triangulated_point);
				next_plps.push_back(next_plp);
				actual_current_plgp = PolyLineGraph2D::plg_point(current_plgp.polyline_id,next_plp);
			} else
				break;
		}
	}

	// Update 3D points
	for(int i=0; i < triangulated_points.size(); i++)
		new_points_data_to_add.push_back(make_pair(triangulated_points[i],PolyLineGraph2D::plg_point(current_plgp.polyline_id,next_plps[i].segment_index,next_plps[i].coords)));

	all_points_compatible = triangulated_points.size() == valid_points.size();

	return triangulated_points.size() > 0;
}

// direction on specified plg, pl compatible with given 3D points? Return amount of compatible points. Doesn't update compatible 3D points
bool compatible_direction_noupdate_vector(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const int current_plg_id, const PolyLineGraph2D::plg_point &current_plgp, const ulong direction, const Mat** all_fundamental_matrices, const vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &valid_points, vector<std::pair<glm::vec3, PolyLineGraph2D::plg_point>> &new_points_data_to_add, bool &all_points_compatible, const int start_interval_index, const int cur_point_index, int end_interval_index, const bool towards_start) {
	new_points_data_to_add = vector<std::pair<glm::vec3, PolyLineGraph2D::plg_point>>();

	const int sz = valid_points.size();
	if(sz == 0)
		return 0;

	bool valid;
	vector<PolyLineGraph2D::polyline::pl_point> next_plps;
	vector<glm::vec3> triangulated_points;

	PolyLineGraph2D::plg_point actual_current_plgp = current_plgp;
	int i= towards_start ? cur_point_index - 1 : cur_point_index + 1;
	while((towards_start && i >= start_interval_index) || (!towards_start && i < end_interval_index)) {
		const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &cur_pt = valid_points[i];

		PolyLineGraph2D::polyline::pl_point next_plp;
		get_plgp_by_epipolar_intersection_from_known_point(sfmd, plgs, current_plg_id, actual_current_plgp, direction, all_fundamental_matrices, cur_pt, next_plp, valid);

		if(!valid)
			break; // fail
		else {
			// check compatibility between new plp and known 3D point observations
			glm::vec3 triangulated_point;
			em_add_new_observation_to_3Dpositions(sfmd,cur_pt,next_plp.coords,current_plg_id,triangulated_point,valid);

			if(valid) {
				triangulated_points.push_back(triangulated_point);
				next_plps.push_back(next_plp);
				actual_current_plgp = PolyLineGraph2D::plg_point(actual_current_plgp.polyline_id,next_plp);
			} else
				break;
		}

		if(towards_start)
			i--;
		else
			i++;

	}

	// Update 3D points
	for(int i=0; i < triangulated_points.size(); i++)
		new_points_data_to_add.push_back(make_pair(triangulated_points[i],PolyLineGraph2D::plg_point(current_plgp.polyline_id,next_plps[i].segment_index,next_plps[i].coords)));

	all_points_compatible = towards_start ? (triangulated_points.size() == cur_point_index) : (triangulated_points.size() == (valid_points.size() - cur_point_index - 1));

	return triangulated_points.size() > 0;
}


// direction on specified plg, pl compatible with given 3D points? If it is compatible with at least 1, return true. Moreover, update compatible 3D points
bool compatible_direction(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const int current_plg_id, const PolyLineGraph2D::plg_point &current_plgp, const ulong direction, const Mat** all_fundamental_matrices, vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &valid_points, bool &all_points_compatible) {
	vector<std::pair<glm::vec3, PolyLineGraph2D::plg_point>> new_points_data_to_add;
	compatible_direction_noupdate(sfmd,plgs,current_plg_id,current_plgp,direction,all_fundamental_matrices,valid_points,new_points_data_to_add,all_points_compatible);

	// Update 3D points
	for(int i=0; i < new_points_data_to_add.size(); i++) {
		std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &to_update = valid_points[i];
		get<0>(to_update) = new_points_data_to_add[i].first;
		get<1>(to_update).push_back(new_points_data_to_add[i].second);
		get<2>(to_update).push_back(current_plg_id);
	}

	return new_points_data_to_add.size() > 0;
}

void find_directions_on_plg_known_3D_point(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const int current_plg_id, const PolyLineGraph2D::plg_point &current_plgp, const Mat** all_fundamental_matrices, vector<ulong> &directions1, const bool direction1_valid, vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &valid_points_direction1, vector<ulong> &directions2, const bool direction2_valid, vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &valid_points_direction2) {
	const PolyLineGraph2D::polyline &pl = plgs[current_plg_id].polylines[current_plgp.polyline_id];
	const ulong start = pl.start;
	const ulong end = pl.end;
	bool all_points_compatible;

	if(direction1_valid) {
		// Try to match direction1 with start
		if(compatible_direction(sfmd, plgs, current_plg_id, current_plgp, start, all_fundamental_matrices,  valid_points_direction1, all_points_compatible)) {
			// direction1 was compatible with start
			directions1[current_plg_id] = start;

			if(direction2_valid) {
				// Check if direction2 is compatible with end
				if(compatible_direction(sfmd, plgs, current_plg_id, current_plgp, end, all_fundamental_matrices,  valid_points_direction2, all_points_compatible))
					directions2[current_plg_id] = end;
			}
		} else if(compatible_direction(sfmd, plgs, current_plg_id, current_plgp, end, all_fundamental_matrices,  valid_points_direction1, all_points_compatible)) {
			// direction1 was compatible with end
			directions1[current_plg_id] = end;

			if(direction2_valid) {
				// Check if direction2 is compatible with start
				if(compatible_direction(sfmd, plgs, current_plg_id, current_plgp, start, all_fundamental_matrices,  valid_points_direction2, all_points_compatible))
					directions2[current_plg_id] = start;
			}
		} else {
			// Couldn't match direction1 in either direction

			if(direction2_valid) {
				if(compatible_direction(sfmd, plgs, current_plg_id, current_plgp, end, all_fundamental_matrices,  valid_points_direction2, all_points_compatible))
					directions2[current_plg_id] = end;
				else if(compatible_direction(sfmd, plgs, current_plg_id, current_plgp, start, all_fundamental_matrices,  valid_points_direction2, all_points_compatible))
					directions2[current_plg_id] = start;
			}
		}
	}
}

void find_directions_on_plg_known_3D_point_no_update(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const int current_plg_id, const PolyLineGraph2D::plg_point &current_plgp, const Mat** all_fundamental_matrices, vector<ulong> &directions1, const bool direction1_valid, const vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &valid_points_direction1, vector<ulong> &directions2, const bool direction2_valid, const vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &valid_points_direction2, vector<std::pair<glm::vec3, PolyLineGraph2D::plg_point>> &new_points_direction1, vector<std::pair<glm::vec3, PolyLineGraph2D::plg_point>> &new_points_direction2) {
	const PolyLineGraph2D::polyline &pl = plgs[current_plg_id].polylines[current_plgp.polyline_id];
	const ulong start = pl.start;
	const ulong end = pl.end;
	bool all_points_compatible;

	if(direction1_valid) {
		// Try to match direction1 with start
		if(compatible_direction_noupdate(sfmd, plgs, current_plg_id, current_plgp, start, all_fundamental_matrices,  valid_points_direction1, new_points_direction1, all_points_compatible)) {
			// direction1 was compatible with start
			directions1[current_plg_id] = start;

			if(direction2_valid) {
				// Check if direction2 is compatible with end
				if(compatible_direction_noupdate(sfmd, plgs, current_plg_id, current_plgp, end, all_fundamental_matrices,  valid_points_direction2, new_points_direction2, all_points_compatible))
					directions2[current_plg_id] = end;
			}
		} else if(compatible_direction_noupdate(sfmd, plgs, current_plg_id, current_plgp, end, all_fundamental_matrices,  valid_points_direction1, new_points_direction1, all_points_compatible)) {
			// direction1 was compatible with end
			directions1[current_plg_id] = end;

			if(direction2_valid) {
				// Check if direction2 is compatible with start
				if(compatible_direction_noupdate(sfmd, plgs, current_plg_id, current_plgp, start, all_fundamental_matrices,  valid_points_direction2, new_points_direction2, all_points_compatible))
					directions2[current_plg_id] = start;
			}
		} else {
			// Couldn't match direction1 in either direction

			if(direction2_valid) {
				if(compatible_direction_noupdate(sfmd, plgs, current_plg_id, current_plgp, end, all_fundamental_matrices,  valid_points_direction2, new_points_direction2, all_points_compatible))
					directions2[current_plg_id] = end;
				else if(compatible_direction_noupdate(sfmd, plgs, current_plg_id, current_plgp, start, all_fundamental_matrices,  valid_points_direction2, new_points_direction2, all_points_compatible))
					directions2[current_plg_id] = start;
			}
		}
	}
}

void find_directions_on_plg_known_3D_point_no_update_vector(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const int current_plg_id, const PolyLineGraph2D::plg_point &current_plgp, const Mat** all_fundamental_matrices, vector<new_3dpoint_plgp_matches> &cur_pts, vector<std::pair<glm::vec3, PolyLineGraph2D::plg_point>> &new_points_direction1, vector<std::pair<glm::vec3, PolyLineGraph2D::plg_point>> &new_points_direction2, ulong &new_direction1, ulong &new_direction2, int start_interval_index, int cur_point_index, int end_interval_index) {
	const PolyLineGraph2D::polyline &pl = plgs[current_plg_id].polylines[current_plgp.polyline_id];
	const ulong start = pl.start;
	const ulong end = pl.end;
	bool all_points_compatible;

	if(cur_point_index > start_interval_index) {
		// Try to match direction1 with start
		//if(compatible_direction_noupdate(sfmd, plgs, current_plg_id, current_plgp, start, all_fundamental_matrices,  valid_points_direction1, new_points_direction1, all_points_compatible)) {
		if(compatible_direction_noupdate_vector(sfmd, plgs, current_plg_id, current_plgp, start, all_fundamental_matrices, cur_pts, new_points_direction1, all_points_compatible, start_interval_index, cur_point_index, end_interval_index, true)) {
			// direction1 was compatible with start
			new_direction1 = start;
			new_direction2 = end;

			if(cur_point_index < end_interval_index) {
				// Check if direction2 is compatible with end
				//if(compatible_direction_noupdate(sfmd, plgs, current_plg_id, current_plgp, end, all_fundamental_matrices,  valid_points_direction2, new_points_direction2, all_points_compatible))
				compatible_direction_noupdate_vector(sfmd, plgs, current_plg_id, current_plgp, end, all_fundamental_matrices, cur_pts, new_points_direction2, all_points_compatible, start_interval_index, cur_point_index, end_interval_index, false);
			}
		//} else if(compatible_direction_noupdate(sfmd, plgs, current_plg_id, current_plgp, end, all_fundamental_matrices,  valid_points_direction1, new_points_direction1, all_points_compatible)) {
		} else if(compatible_direction_noupdate_vector(sfmd, plgs, current_plg_id, current_plgp, end, all_fundamental_matrices, cur_pts, new_points_direction1, all_points_compatible, start_interval_index, cur_point_index, end_interval_index, true)) {
			// direction1 was compatible with end
			new_direction1 = end;
			new_direction2 = start;

			if(cur_point_index < end_interval_index) {
				// Check if direction2 is compatible with start
				//if(compatible_direction_noupdate(sfmd, plgs, current_plg_id, current_plgp, start, all_fundamental_matrices,  valid_points_direction2, new_points_direction2, all_points_compatible))
				compatible_direction_noupdate_vector(sfmd, plgs, current_plg_id, current_plgp, start, all_fundamental_matrices, cur_pts, new_points_direction2, all_points_compatible, start_interval_index, cur_point_index, end_interval_index, false);

			}
		} else {
			// Couldn't match direction1 in either direction

			if(cur_point_index < end_interval_index) {
				//if(compatible_direction_noupdate(sfmd, plgs, current_plg_id, current_plgp, end, all_fundamental_matrices,  valid_points_direction2, new_points_direction2, all_points_compatible))
				if(compatible_direction_noupdate_vector(sfmd, plgs, current_plg_id, current_plgp, end, all_fundamental_matrices, cur_pts, new_points_direction2, all_points_compatible, start_interval_index, cur_point_index, end_interval_index, false)) {
					new_direction2 = end;
					new_direction1 = start;
				// else if(compatible_direction_noupdate(sfmd, plgs, current_plg_id, current_plgp, start, all_fundamental_matrices,  valid_points_direction2, new_points_direction2, all_points_compatible))
			}else if(compatible_direction_noupdate_vector(sfmd, plgs, current_plg_id, current_plgp, start, all_fundamental_matrices, cur_pts, new_points_direction2, all_points_compatible, start_interval_index, cur_point_index, end_interval_index, false)) {
					new_direction2 = start;
					new_direction1 = end;
			}
			}
		}
	}
}

void find_directions_all_views(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs, const std::tuple<glm::vec3,vector<PolyLineGraph2D::plg_point>,vector<int>> &matches, vector<ulong> &direction1, bool &direction1_valid, vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &valid_points_direction1, vector<ulong> &direction2, bool &direction2_valid, vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &valid_points_direction2, bool &fail_due_to_parallel_epipolar) {
	int i;

	find_directions_3view_firstlast(sfmd, all_fundamental_matrices, plgs, matches, direction1, direction1_valid, valid_points_direction1, direction2, direction2_valid, valid_points_direction2, fail_due_to_parallel_epipolar);

	if(direction1_valid) {
		const vector<PolyLineGraph2D::plg_point> &matches_data = get<1>(matches);
		const vector<int> &matched_ids = get<2>(matches);
		const int amount_of_matches = matched_ids.size();

		for(i=1; i < amount_of_matches / 2; i++)
			find_directions_on_plg_known_3D_point(sfmd,plgs,matched_ids[i], matches_data[i], all_fundamental_matrices,direction1,direction1_valid,valid_points_direction1,direction2,direction2_valid,valid_points_direction2);

		for(i=amount_of_matches / 2 + 1; i < amount_of_matches-1; i++)
			find_directions_on_plg_known_3D_point(sfmd,plgs,matched_ids[i], matches_data[i],all_fundamental_matrices,direction1,direction1_valid,valid_points_direction1,direction2,direction2_valid,valid_points_direction2);
	}
}

/*
*
 * Follow given directions on all plgs until valid.

void follow_direction_second(const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const vector<ulong> &directions, const Mat** all_fundamental_matrices, vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &valid_points, bool &fail_due_to_parallel_epipolar) {
	std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> new_point_data;
	// cout << "Following direction\n";
	while(compatible_second(sfmd, plgs, directions, all_fundamental_matrices, valid_points[valid_points.size()-1], new_point_data, fail_due_to_parallel_epipolar))
		valid_points.push_back(std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>(new_point_data));
}*/

pair<vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>>,vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>>> follow_plgs_from_match2(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs, const std::tuple<glm::vec3,vector<PolyLineGraph2D::plg_point>,vector<int>> &matches, bool &valid) {
	const vector<PolyLineGraph2D::plg_point> &matches_data = get<1>(matches);
	const vector<int> &matched_ids = get<2>(matches);
	const int amount_of_matches = matched_ids.size();
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> points_dir1, points_dir2;
	bool fail_due_to_parallel_epipolar;
	valid = false;

	if(amount_of_matches < 3)
		return make_pair(points_dir1,points_dir2);

	int selected_indexes[3];

	selected_indexes[0] = 0;
	selected_indexes[1] = amount_of_matches/2;
	selected_indexes[2] = amount_of_matches-1;

	vector<int> selected_2d_reprojections_ids;
	selected_2d_reprojections_ids.push_back(matched_ids[selected_indexes[0]]);
	selected_2d_reprojections_ids.push_back(matched_ids[selected_indexes[1]]);
	selected_2d_reprojections_ids.push_back(matched_ids[selected_indexes[2]]);

	const tuple<PolyLineGraph2DHMapImpl,PolyLineGraph2DHMapImpl,PolyLineGraph2DHMapImpl> selected_plgs = make_tuple(plgs[selected_2d_reprojections_ids[0]],plgs[selected_2d_reprojections_ids[1]],plgs[selected_2d_reprojections_ids[2]]);
	const tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point> current_plgps = make_tuple(matches_data[selected_indexes[0]],matches_data[selected_indexes[1]],matches_data[selected_indexes[2]]);
	tuple<ulong,ulong,ulong> direction1,direction2;
	bool direction1_valid,direction2_valid;
	vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> valid_points_direction1,valid_points_direction2;

	const Mat &fab = all_fundamental_matrices[selected_2d_reprojections_ids[0]][selected_2d_reprojections_ids[1]];
	const Mat &fac = all_fundamental_matrices[selected_2d_reprojections_ids[0]][selected_2d_reprojections_ids[2]];

	find_directions(sfmd, selected_2d_reprojections_ids, selected_plgs, current_plgps, fab, fac, direction1, direction1_valid, valid_points_direction1, direction2, direction2_valid, valid_points_direction2,fail_due_to_parallel_epipolar);

	if(direction1_valid) {
		for(const auto &p : valid_points_direction1) {
			const tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point> &p_plgps = p.first;
			const std::tuple<glm::vec3, vector<glm::vec2>, vector<int>> &p3d = p.second;

			vector<PolyLineGraph2D::plg_point> selected_plgps;
			selected_plgps.push_back(get<0>(p_plgps));
			selected_plgps.push_back(get<1>(p_plgps));
			selected_plgps.push_back(get<2>(p_plgps));

			points_dir1.push_back(make_tuple(get<0>(p3d),selected_plgps,get<2>(p3d)));
		}

		vector<ulong> directions1(plgs.size());
		directions1[get<2>(valid_points_direction1[0].second)[0]] = get<0>(direction1);
		directions1[get<2>(valid_points_direction1[0].second)[1]] = get<1>(direction1);
		directions1[get<2>(valid_points_direction1[0].second)[2]] = get<2>(direction1);

		follow_direction(sfmd, plgs, directions1, all_fundamental_matrices, points_dir1, fail_due_to_parallel_epipolar);
	}

	if(direction2_valid) {
		for(const auto &p : valid_points_direction2) {
			const tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point> &p_plgps = p.first;
			const std::tuple<glm::vec3, vector<glm::vec2>, vector<int>> &p3d = p.second;

			vector<PolyLineGraph2D::plg_point> selected_plgps;
			selected_plgps.push_back(get<0>(p_plgps));
			selected_plgps.push_back(get<1>(p_plgps));
			selected_plgps.push_back(get<2>(p_plgps));

			points_dir2.push_back(make_tuple(get<0>(p3d),selected_plgps,get<2>(p3d)));
		}

		vector<ulong> directions2(plgs.size());
		directions2[get<2>(valid_points_direction2[0].second)[0]] = get<0>(direction2);
		directions2[get<2>(valid_points_direction2[0].second)[1]] = get<1>(direction2);
		directions2[get<2>(valid_points_direction2[0].second)[2]] = get<2>(direction2);

		follow_direction(sfmd, plgs, directions2, all_fundamental_matrices, points_dir2, fail_due_to_parallel_epipolar);
	}

	return make_pair(points_dir1,points_dir2);
}

pair<vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>>,vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>>> follow_plgs_from_match3(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs, const std::tuple<glm::vec3,vector<PolyLineGraph2D::plg_point>,vector<int>> &matches, bool &valid) {
	const vector<PolyLineGraph2D::plg_point> &matches_data = get<1>(matches);
	const vector<int> &matched_ids = get<2>(matches);
	const int amount_of_matches = matched_ids.size();
	bool fail_due_to_parallel_epipolar;
	valid = false;

	if(amount_of_matches < 3)
		return make_pair(vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>>(),vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>>());

	bool direction1_valid,direction2_valid;

	vector<ulong> direction1,direction2;
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> valid_points_direction1,valid_points_direction2;

	find_directions_all_views(sfmd, all_fundamental_matrices, plgs, matches, direction1, direction1_valid, valid_points_direction1, direction2, direction2_valid, valid_points_direction2, fail_due_to_parallel_epipolar);

	// cout << "Direction 1 valid: " << direction1_valid<< endl;
	// cout << "Direction 2 valid: " << direction2_valid<< endl;

	// cout << (direction2_valid || direction1_valid ? ">>>>>>> Polyline follow : SUCCESS\n" : ">>>>>>> Polyline follow : FAIL\n");

	if(direction1_valid)
		follow_direction(sfmd, plgs, direction1, all_fundamental_matrices, valid_points_direction1, fail_due_to_parallel_epipolar);


	if(direction2_valid)
		follow_direction(sfmd, plgs, direction2, all_fundamental_matrices, valid_points_direction2, fail_due_to_parallel_epipolar);


	return make_pair(valid_points_direction1,valid_points_direction2);
}

vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> follow_plgs_from_match(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs, const std::tuple<glm::vec3,vector<PolyLineGraph2D::plg_point>,vector<int>> &matches, bool &valid) {
	const vector<PolyLineGraph2D::plg_point> &matches_data = get<1>(matches);
	const vector<int> &matched_ids = get<2>(matches);
	const int amount_of_matches = matched_ids.size();
	vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> valid_points;
	bool fail_due_to_parallel_epipolar;
	valid = false;

	if(amount_of_matches < 3) {
		// cout << "Less than 3 initial matches, cannot follow plgpoint!\n";
		return valid_points;
	}

	int selected_indexes[3];

	selected_indexes[0] = 0;
	selected_indexes[1] = amount_of_matches/2;
	selected_indexes[2] = amount_of_matches-1;

	vector<int> selected_2d_reprojections_ids;
	selected_2d_reprojections_ids.push_back(matched_ids[selected_indexes[0]]);
	selected_2d_reprojections_ids.push_back(matched_ids[selected_indexes[1]]);
	selected_2d_reprojections_ids.push_back(matched_ids[selected_indexes[2]]);

	const tuple<PolyLineGraph2DHMapImpl,PolyLineGraph2DHMapImpl,PolyLineGraph2DHMapImpl> selected_plgs = make_tuple(plgs[selected_2d_reprojections_ids[0]],plgs[selected_2d_reprojections_ids[1]],plgs[selected_2d_reprojections_ids[2]]);
	const tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point> current_plgps = make_tuple(matches_data[selected_indexes[0]],matches_data[selected_indexes[1]],matches_data[selected_indexes[2]]);
	tuple<ulong,ulong,ulong> direction1,direction2;
	bool direction1_valid,direction2_valid;
	vector<pair<tuple<PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point,PolyLineGraph2D::plg_point>,std::tuple<glm::vec3, vector<glm::vec2>, vector<int>>>> valid_points_direction1,valid_points_direction2;

	const Mat &fab = all_fundamental_matrices[selected_2d_reprojections_ids[0]][selected_2d_reprojections_ids[1]];
	const Mat &fac = all_fundamental_matrices[selected_2d_reprojections_ids[0]][selected_2d_reprojections_ids[2]];

	find_directions(sfmd, selected_2d_reprojections_ids, selected_plgs, current_plgps, fab, fac, direction1, direction1_valid, valid_points_direction1, direction2, direction2_valid, valid_points_direction2, fail_due_to_parallel_epipolar);

	if(direction1_valid)
		for(const auto &vp : valid_points_direction1)
			valid_points.push_back(vp);

	if(direction1_valid)
		for(const auto &vp : valid_points_direction2)
			valid_points.push_back(vp);

	return valid_points;
}




void follow_plgs_from_match4(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs, const new_3dpoint_plgp_matches &matches, vector<ulong> &directions1, bool &direction1_valid, vector<new_3dpoint_plgp_matches> &valid_points_direction1, vector<ulong> &directions2, bool &direction2_valid, vector<new_3dpoint_plgp_matches> &valid_points_direction2) {
	const vector<PolyLineGraph2D::plg_point> &matches_data = get<1>(matches);
	const vector<int> &matched_ids = get<2>(matches);
	const int amount_of_matches = matched_ids.size();
	bool fail_due_to_parallel_epipolar;

	direction1_valid=false;
	direction2_valid=false;

	if(amount_of_matches < 3)
		return;

	find_directions_all_views(sfmd, all_fundamental_matrices, plgs, matches, directions1, direction1_valid, valid_points_direction1, directions2, direction2_valid, valid_points_direction2, fail_due_to_parallel_epipolar);

	if(direction1_valid)
		follow_direction(sfmd, plgs, directions1, all_fundamental_matrices, valid_points_direction1, fail_due_to_parallel_epipolar);


	if(direction2_valid)
		follow_direction(sfmd, plgs, directions2, all_fundamental_matrices, valid_points_direction2, fail_due_to_parallel_epipolar);

}


/**
 * A new given 3D plg point is compatible if it is possible to discover at least 2 new points in a direction by following plgs
 */
bool compatible_new_plg_point(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs, const new_3dpoint_plgp_matches &matches, vector<ulong> &directions1, bool &direction1_valid, vector<new_3dpoint_plgp_matches> &valid_points_direction1, vector<ulong> &directions2, bool &direction2_valid, vector<new_3dpoint_plgp_matches> &valid_points_direction2) {

	follow_plgs_from_match4(sfmd, all_fundamental_matrices, plgs, matches, directions1, direction1_valid, valid_points_direction1, directions2, direction2_valid, valid_points_direction2);

	// previous
	if(direction1_valid && valid_points_direction1.size() >= 2)
		return true;
	if(direction2_valid && valid_points_direction2.size() >= 2)
		return true;

	return false;
}

bool add_view_to_3dpoint_and_sides_plgp_matches(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs, new_3dpoint_and_sides_plgp_matches &cur_pts, const int current_plg_id, const PolyLineGraph2D::plg_point &current_plgp) {
	// check if central point can be added
	glm::vec3 new_central_point_coords;
	if(!compatible_new_observation_to_3Dpositions(sfmd, get<1>(cur_pts), current_plgp, current_plg_id,new_central_point_coords))
		return false;

	vector<std::pair<glm::vec3, PolyLineGraph2D::plg_point>> new_points_direction1;
	vector<std::pair<glm::vec3, PolyLineGraph2D::plg_point>> new_points_direction2;
	vector<ulong> &directions1 = get<0>(cur_pts).second;
	vector<ulong> &directions2 = get<2>(cur_pts).second;

	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &valid_points_direction1 = get<0>(cur_pts).first;
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &valid_points_direction2 = get<2>(cur_pts).first;

	bool direction1_valid = valid_points_direction1.size() > 0;
	bool direction2_valid = valid_points_direction2.size() > 0;


	find_directions_on_plg_known_3D_point_no_update(sfmd, plgs, current_plg_id, current_plgp, all_fundamental_matrices, directions1, direction1_valid, valid_points_direction1, directions2, direction2_valid, valid_points_direction2,new_points_direction1,new_points_direction2);

#if defined(SWITCH_PLG_MATCHING_ADDPOINT_BOTHDIR_ONE)
	if(direction1_valid && new_points_direction1.size() == 0)
		return false; // fail

	if(direction2_valid && new_points_direction2.size() == 0)
		return false; // fail
#endif

	// Success! update 3D points

	update_new_3dpoint_plgp_matches(get<1>(cur_pts),current_plg_id,current_plgp,new_central_point_coords);

	for(int i=0; i < new_points_direction1.size(); i++)
		update_new_3dpoint_plgp_matches(valid_points_direction1[i], current_plg_id,new_points_direction1[i].second,	new_points_direction1[i].first);

	for(int i=0; i < new_points_direction2.size(); i++)
		update_new_3dpoint_plgp_matches(valid_points_direction2[i],current_plg_id,new_points_direction2[i].second,new_points_direction2[i].first);

	bool fail_due_to_parallel_epipolar;

	if(direction1_valid && new_points_direction1.size() == valid_points_direction1.size())
		follow_direction(sfmd, plgs, directions1, all_fundamental_matrices, valid_points_direction1, fail_due_to_parallel_epipolar);

	if(direction2_valid && new_points_direction2.size() == valid_points_direction2.size())
		follow_direction(sfmd, plgs, directions2, all_fundamental_matrices, valid_points_direction2, fail_due_to_parallel_epipolar);


	return true;
}

// Try to add current PLGP to the current matched points by Polyline following
// Does not go outside bounds [start_interval_index,end_interval_index)
// cur_plgp is first matched with cur_p3ds[start_interval_index]
// then tries to follow polyine and match in both directions
// in the direction of the start of cur_p3ds, doesn't go over start_interval_index (unless start_interval_index is zero, in which case it tries to find new points in that direction as normal)
// Return pair<amount_of_matches_towards_start,amount_of_matches_towards_end>
pair<int,int> add_view_to_3dpoint_and_sides_plgp_matches_vector(const SfMData &sfmd, const Mat** all_fundamental_matrices, const vector<PolyLineGraph2DHMapImpl> &plgs, vector<new_3dpoint_plgp_matches> &cur_pts, vector<ulong> &start_dirs, vector<ulong> &end_dirs, const int current_plg_id, const PolyLineGraph2D::plg_point &current_plgp, int start_interval_index, int cur_point_index, int end_interval_index, bool &success) {
	success = false;

	// check if central point can be added
	glm::vec3 new_central_point_coords;
	if(!compatible_new_observation_to_3Dpositions(sfmd, cur_pts[cur_point_index], current_plgp, current_plg_id,new_central_point_coords))
		return make_pair(0,0);

	ulong new_direction1,new_direction2;

	int amount_of_matches_towards_start = 0;
	int amount_of_matches_towards_end =0;

	vector<std::pair<glm::vec3, PolyLineGraph2D::plg_point>> new_points_direction1;
	vector<std::pair<glm::vec3, PolyLineGraph2D::plg_point>> new_points_direction2;

	find_directions_on_plg_known_3D_point_no_update_vector(sfmd, plgs, current_plg_id, current_plgp, all_fundamental_matrices, cur_pts, new_points_direction1, new_points_direction2, new_direction1, new_direction2, start_interval_index, cur_point_index, end_interval_index);

#if defined(SWITCH_PLG_MATCHING_ADDPOINT_BOTHDIR_ONE)
	if(cur_point_index > 0 && new_points_direction1.size() == 0)
		return make_pair(0,0); // fail

	if(cur_point_index < cur_pts.size()-1 && new_points_direction2.size() == 0)
		return make_pair(0,0); // fail
#else
	return make_pair(0,0); // fail
#endif

	// SUCCESS

	success = true;

	amount_of_matches_towards_start = new_points_direction1.size();
	amount_of_matches_towards_end = new_points_direction2.size();

	update_new_3dpoint_plgp_matches(cur_pts[cur_point_index],current_plg_id,current_plgp,new_central_point_coords);

	for(int i=0; i < new_points_direction1.size(); i++)
		update_new_3dpoint_plgp_matches(cur_pts[cur_point_index-1-i], current_plg_id, new_points_direction1[i].second, new_points_direction1[i].first);

	for(int i=0; i < new_points_direction2.size(); i++)
		update_new_3dpoint_plgp_matches(cur_pts[cur_point_index+1+i],current_plg_id,new_points_direction2[i].second,new_points_direction2[i].first);

	bool fail_due_to_parallel_epipolar;
	int starting_sz;

	int amount_of_new_matches_towards_start = 0;

	if(new_points_direction1.size() > 0 && new_points_direction1.size() == cur_point_index) {
		// expand by follow towards start
		starting_sz = cur_pts.size();
		start_dirs[current_plg_id] = new_direction1;
		follow_direction_vector_start(sfmd, plgs, start_dirs, all_fundamental_matrices, cur_pts, fail_due_to_parallel_epipolar);
		amount_of_new_matches_towards_start = (cur_pts.size() - starting_sz);
		amount_of_matches_towards_start += amount_of_new_matches_towards_start;
		cur_point_index += amount_of_new_matches_towards_start;
	}

	if(new_points_direction2.size() > 0 && new_points_direction2.size() == (cur_pts.size() - cur_point_index - 1)) {
		// expand by follow towards end
		starting_sz = cur_pts.size();
		end_dirs[current_plg_id] = new_direction2;
		follow_direction_vector_end(sfmd, plgs, end_dirs, all_fundamental_matrices, cur_pts, fail_due_to_parallel_epipolar);
		amount_of_matches_towards_end += (cur_pts.size() - starting_sz);
	}

	return make_pair(amount_of_matches_towards_start,amount_of_matches_towards_end);
}



