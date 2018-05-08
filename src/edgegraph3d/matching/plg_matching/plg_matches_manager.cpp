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


#include "plg_matches_manager.hpp"

#include <fstream>

#include <iostream>
#include <map>
#include <utility>

#include "geometric_utilities.hpp"

PLGMatchesManager::PLGMatchesManager(vector<PolyLineGraph2DHMapImpl> &plgs, PolyLineGraph3DHMapImpl &plg3d) : plgs(plgs), plg3d(plg3d) {
	omp_init_lock(&writelock); // Initialize write lock
	for(int plg_id=0; plg_id < plgs.size(); plg_id++)
		matched_polyline_intervals.push_back(vector<sorted_pl_intervals_set>(plgs[plg_id].get_polylines_amount()));
}

PLGMatchesManager::~PLGMatchesManager() {}

PLGMatchesManager::PLGMatchesManager() {}

PLGMatchesManager::PLGMatchesManager(vector<PolyLineGraph2DHMapImpl> &plgs, PolyLineGraph3DHMapImpl &plg3d, vector<PolyLineGraph2D::plgs_point> &point_matches, plgsp_to_p3d_id_map &point_matches_map, pointmap3dtoplgpstype &p3d_to_matches_map, vector<vector<sorted_pl_intervals_set>> &matched_polyline_intervals) :
	plgs(plgs), plg3d(plg3d), point_matches(point_matches), point_matches_map(point_matches_map), matched_polyline_intervals(matched_polyline_intervals) {}

bool PLGMatchesManager::is_matched(const int plg_id,const ulong polyline_id, const PolyLineGraph2D::polyline::pl_point &plp, PolyLineGraph2D::polyline::pl_interval &pli_containing_point) {
	//return false;
	sorted_pl_intervals_set& intervals = matched_polyline_intervals[plg_id][polyline_id];

	if(intervals.size() > 0) {
		const PolyLineGraph2D::polyline &pl = plgs[plg_id].polylines[polyline_id];
		sorted_pl_intervals_set::iterator it = intervals.lower_bound(PolyLineGraph2D::polyline::pl_interval(plp,PolyLineGraph2D::polyline::pl_point()));

		if(it == intervals.end() || it->start.segment_index > plp.segment_index)
			it--;

		// Now it is on an element with a starting index smaller or equal to plp.segment_index
		bool reached_begin=false;
		bool found = false;
		for(;it->start.segment_index == plp.segment_index && !reached_begin; it--) {
			if(it == intervals.begin())
				reached_begin = true;
			if(pl.interval_contains_plp(*it,plp)) {
				found = true;
				pli_containing_point = *it;
				break;
			}
		}
		if(!found && !reached_begin && it->end.segment_index <= plp.segment_index)
			if(pl.interval_contains_plp(*it,plp)) {
				found = true;
				pli_containing_point = *it;
			}
		return found;
	} else
		return false;
}

bool PLGMatchesManager::is_matched(const int plg_id, const PolyLineGraph2D::plg_point &plgp, PolyLineGraph2D::polyline::pl_interval &pli_containing_point) {
	return is_matched(plg_id, plgp.polyline_id, plgp.plp,pli_containing_point);
}

bool PLGMatchesManager::is_matched(const PolyLineGraph2D::plgs_point &plgsp, PolyLineGraph2D::polyline::pl_interval &pli_containing_point) {
	return is_matched(plgsp.plg_id,plgsp.plgp,pli_containing_point);
}

PolyLineGraph3DHMapImpl& PLGMatchesManager::get_plg3d() {
	return plg3d;
}

void PLGMatchesManager::add_matched_2dsegment(const int plg_id, const ulong pl_id, const PolyLineGraph2D::polyline::pl_point &pl_a, const PolyLineGraph2D::polyline::pl_point &pl_b) {
	if(pl_a.segment_index < pl_b.segment_index)
		matched_polyline_intervals[plg_id][pl_id].insert(PolyLineGraph2D::polyline::pl_interval(pl_a,pl_b));
	else if (pl_a.segment_index > pl_b.segment_index)
		matched_polyline_intervals[plg_id][pl_id].insert(PolyLineGraph2D::polyline::pl_interval(pl_b,pl_a));
	else if (is_ordered_2dlinepoints(plgs[plg_id].polylines[pl_id].polyline_coords[pl_a.segment_index],pl_a.coords,pl_b.coords))
		matched_polyline_intervals[plg_id][pl_id].insert(PolyLineGraph2D::polyline::pl_interval(pl_a,pl_b));
	else
		matched_polyline_intervals[plg_id][pl_id].insert(PolyLineGraph2D::polyline::pl_interval(pl_b,pl_a));
}

void PLGMatchesManager::add_matched_3dsegment(const new_3dpoint_plgp_matches &p1, const new_3dpoint_plgp_matches &p2) {
	pair<ulong,ulong> node_ids;
	plg3d.add_direct_connection(get<0>(p1), get<0>(p2), node_ids);

	// adding 2D observations
	plg3d.set_observations(node_ids.first, convert_vecplgp_to_vec2(get<1>(p1)), get<2>(p1));
	plg3d.set_observations(node_ids.second, convert_vecplgp_to_vec2(get<1>(p2)), get<2>(p2));

	p3d_to_matches_map[get<0>(p1)]=p1;
	p3d_to_matches_map[get<0>(p2)]=p2;

	vector<PolyLineGraph2D::plg_point> plgps1(plgs.size());
	vector<bool> seen1(plgs.size(),false);
	vector<PolyLineGraph2D::plg_point> plgps2(plgs.size());
	vector<bool> seen2(plgs.size(),false);


	for(int i=0; i < get<1>(p1).size(); i++) {
		seen1[get<2>(p1)[i]] = true;
		plgps1[get<2>(p1)[i]]=get<1>(p1)[i];

		if(!plgs[get<2>(p1)[i]].is_valid_polyline(get<1>(p1)[i].polyline_id))
			cout << "Invalid polyline selected!\n";
	}

	for(int i=0; i < get<1>(p2).size(); i++) {
		seen2[get<2>(p2)[i]] = true;
		plgps2[get<2>(p2)[i]]=get<1>(p2)[i];

		if(!plgs[get<2>(p2)[i]].is_valid_polyline(get<1>(p2)[i].polyline_id))
			cout << "Invalid polyline selected!\n";
	}

	for(int plg_id=0; plg_id < plgs.size(); plg_id++)
		if(seen1[plg_id] && seen2[plg_id]) {
			if(plgps1[plg_id].polyline_id != plgps2[plg_id].polyline_id) {
				// Then either the first point is extreme of its polyline connected to the next polyline or there is a problem
				const PolyLineGraph2D::polyline &pl = plgs[plg_id].polylines[plgps1[plg_id].polyline_id];
				ulong node_id;
				bool valid = false;
				bool is_extreme = false;
				if(pl.is_start(plgps1[plg_id].plp)) {
					node_id = pl.start;
					is_extreme = true;
				}
				if(!is_extreme && pl.is_end(plgps1[plg_id].plp)) {
					node_id = pl.end;
					is_extreme = true;
				}
				if(is_extreme) {
					const PolyLineGraph2D::polyline &next_pl = plgs[plg_id].polylines[plgps2[plg_id].polyline_id];
					PolyLineGraph2D::polyline::pl_point extreme_second_pl = next_pl.get_extreme_plp(node_id,valid);
					if(valid) // If the two polylines are actually connected at node_id
						add_matched_2dsegment(plg_id,plgps2[plg_id].polyline_id,extreme_second_pl,plgps2[plg_id].plp);
				}

/*				// Commented because of expand all views
 * 				else
					std::invalid_argument("add_matched_3dsegment - trying to add pl_interval on different polylines");*/

			} else
				add_matched_2dsegment(plg_id,plgps1[plg_id].polyline_id,plgps1[plg_id].plp,plgps2[plg_id].plp);
		}
}

void PLGMatchesManager::add_matched_3dpolyline(const vector<new_3dpoint_plgp_matches> &pl) {
	omp_set_lock(&writelock);
	for(int i=1; i < pl.size(); i++)
		add_matched_3dsegment(pl[i-1],pl[i]);
	omp_unset_lock(&writelock);
}

void serialize_plgmm(const PLGMatchesManager &plgmm, const string plgmm_path) {
	std::ofstream ofs(plgmm_path);
    boost::archive::text_oarchive oa(ofs);
    oa << plgmm;
}

PLGMatchesManager deserialize_plgmm(const string plgmm_path) {
	PLGMatchesManager plgmm;
    std::ifstream ifs(plgmm_path);
    boost::archive::text_iarchive ia(ifs);
    ia >> plgmm;
    return plgmm;
}

