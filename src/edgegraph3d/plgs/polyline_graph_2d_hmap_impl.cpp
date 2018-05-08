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

#include <iostream>
#include <map>
#include <set>
#include <utility>

#include "polyline_graph_2d_hmap_impl.hpp"
#include "datatypes.hpp"
#include "geometric_utilities.hpp"
#include "glm.hpp"


PolyLineGraph2DHMapImpl::PolyLineGraph2DHMapImpl() {
}

PolyLineGraph2DHMapImpl::PolyLineGraph2DHMapImpl(const vector<polyline> &polylines, const vector<vector<ulong>> &connections, const vector<bool> &visited_nodes, const ulong &real_nodes_amount, const ulong &nodes_amount, const ulong &next_node_id, const vector<glm::vec2> &nodes_coords, const pointmaptype &point_map) :
	PolyLineGraph2D(polylines, connections, visited_nodes, real_nodes_amount, nodes_amount, next_node_id, nodes_coords), point_map(point_map) {}

PolyLineGraph2DHMapImpl::~PolyLineGraph2DHMapImpl() {
	//delete point_map;
}

ulong PolyLineGraph2DHMapImpl::get_node_id(const glm::vec2 &p_coords) {
	pointmaptype::const_iterator it = point_map.find(p_coords);
	ulong node_id;

	if(it != point_map.end() && !is_valid_node(it->second)) {
		invalidate_node(it->second);
		it = point_map.end();
	}

	if(it == point_map.end()) {
		// Create new node, since the coordinates are new
		//node_id = point_map.size();
		node_id = get_next_node_id();
		point_map[p_coords] = node_id;
		connections.push_back(vector<ulong>());
		PolyLineGraph2D::add_node_coords(p_coords);
		increment_real_nodes_amount();
	} else
		node_id = it->second;

	return node_id;
}

void PolyLineGraph2DHMapImpl::internal_add_polyline(const polyline &pl) {
	if(!is_duplicate(pl)) {
		ulong pl_id = polylines.size();
		polylines.push_back(pl);
		connections[pl.start].push_back(pl_id);
		if(pl.start != pl.end)
			connections[pl.end].push_back(pl_id);
	}
}

vector<glm::vec2> filter_polyline(const vector<glm::vec2> &polyline_coords) {
	vector<glm::vec2> new_plc;
	if(polyline_coords[0] == polyline_coords[polyline_coords.size()-1] && polyline_coords.size() == 4 && squared_2d_distance(polyline_coords[1],polyline_coords[2]) <= 4) {
		new_plc.push_back(polyline_coords[0]);
		new_plc.push_back(middle_point(polyline_coords[1],polyline_coords[2]));
	} else
		new_plc = polyline_coords;
	return new_plc;
}

/**
 * Note: check if identical polyline is already in
 */
void PolyLineGraph2DHMapImpl::add_polyline(const vector<glm::vec2> &polyline_coords) {
	const vector<glm::vec2> polyline_coords2 = filter_polyline(polyline_coords);
	ulong start = get_node_id(polyline_coords2[0]);
	ulong end = get_node_id(polyline_coords2[polyline_coords2.size()-1]);
	polyline pl(start,end,polyline_coords2);
	internal_add_polyline(pl);
}

ulong PolyLineGraph2DHMapImpl::split_polyline(const PolyLineGraph2D::plg_point &plgp) {
	ulong split_node_id = get_node_id(plgp.plp.coords);
	pair<vector<glm::vec2>,vector<glm::vec2>> pls = polylines[plgp.polyline_id].split(plgp.plp);

	// NOTE: putting remove_polyline before the add_polyline s seems to be bugged
	remove_polyline(plgp.polyline_id);

	add_polyline(pls.first);
	add_polyline(pls.second);

	return split_node_id;
}

bool PolyLineGraph2DHMapImpl::is_duplicate(const polyline &pl) {
	const vector<ulong> &s_connections = connections[pl.start];
	const vector<ulong> &e_connections = connections[pl.end];

	const vector<ulong> &smallest_connections = s_connections.size() < e_connections.size() ? s_connections : e_connections;

	for(const auto pl_id : smallest_connections)
		if(polylines[pl_id] == pl)
			return true;

	return false;
}

void PolyLineGraph2DHMapImpl::add_direct_connection(const ulong start, const ulong end) {
	vector<glm::vec2> connection;
	connection.push_back(nodes_coords[start]);
	connection.push_back(nodes_coords[end]);
	polyline pl(start,end,connection);
	internal_add_polyline(pl);
}



void PolyLineGraph2DHMapImpl::connect_close_extremes() {
	vector<pair<ulong,ulong>> extreme_pairs = find_closest_pairs_with_max_dist(get_extreme_nodes_ids_and_coords(),DIRECT_CONNECTION_EXTREMES_MAXDIST);

	pair<vector<ulong>,vector<set<ulong>>> components = compute_components();
	vector<ulong> node_components_ids = components.first;
	vector<set<ulong>> components_nodes = components.second;

	ulong component_to_change;
	ulong new_id;

	for(const auto &pp : extreme_pairs)
		if(node_components_ids[std::get<0>(pp)] != node_components_ids[std::get<1>(pp)]) {
			if(intersect_polylines(glm::vec4(nodes_coords[std::get<0>(pp)],nodes_coords[std::get<1>(pp)])).size() == 0) {
				add_direct_connection(std::get<0>(pp),std::get<1>(pp));
				if(components_nodes[node_components_ids[std::get<0>(pp)]].size() < components_nodes[node_components_ids[std::get<1>(pp)]].size()) {
					new_id = node_components_ids[std::get<1>(pp)];
					component_to_change = node_components_ids[std::get<0>(pp)];
				} else {
					new_id = node_components_ids[std::get<0>(pp)];
					component_to_change = node_components_ids[std::get<1>(pp)];
				}
				set<ulong> &to_change = components_nodes[component_to_change];
				for(std::set<ulong>::iterator it = to_change.begin(); it != to_change.end(); it++) {
					node_components_ids[*it] = new_id;
				}
			}
		}
}

void PolyLineGraph2DHMapImpl::invalidate_node(ulong node_id) {
	point_map.erase(nodes_coords[node_id]);
	PolyLineGraph2D::invalidate_node(node_id);
}

void PolyLineGraph2DHMapImpl::remove_2connection_nodes() {
	for(ulong node_id = 0; node_id < connections.size(); node_id++)
		if(connections[node_id].size() == 2) {
			// cout << "Found 2-connect node " << node_id << " with coords " << nodes_coords[node_id] << " : ";
			const ulong p_id_1 = connections[node_id][0];
			polyline &p1 = polylines[p_id_1];
			const ulong p_id_2 = connections[node_id][1];
			polyline &p2 = polylines[p_id_2];

			const ulong other_end_1 = p1.get_other_end(node_id);
			const ulong other_end_2 = p2.get_other_end(node_id);

			if(vec_glm_vec2_equal(p1.polyline_coords,p2.polyline_coords) || vec_glm_vec2_equal_inv(p1.polyline_coords,p2.polyline_coords)) {
				remove_polyline(p_id_2);
				continue;
			}

			if(other_end_1 != node_id && other_end_2 != node_id) {
				polyline p3 = PolyLineGraph2D::polyline::merge_polylines(p1,p2);

				internal_add_polyline(p3);
				remove_polyline(p_id_1);
				remove_polyline(p_id_2);
				invalidate_node(node_id);
			}
		}
}

void PolyLineGraph2DHMapImpl::remove_degenerate_loops() {
	const ulong polylines_amount = polylines.size();
	for(ulong i=0; i < polylines_amount; i++)
		if(is_valid_polyline(i)) {
			const PolyLineGraph2D::polyline &p = polylines[i];
			if(p.start == p.end || p.polyline_coords[0] == p.polyline_coords[p.polyline_coords.size()-1]) {
				if(p.polyline_coords.size() < 5)
					remove_polyline(i);
			}
		}
}

void PolyLineGraph2DHMapImpl::remove_invalid_polylines() {
	const ulong polylines_amount = polylines.size();
	for(ulong i=0; i < polylines_amount; i++)
		if(!is_valid_polyline(i))
			remove_polyline(i);

}

#define SHORT_POLYLINES_MAXLENGTH 2.0

void PolyLineGraph2DHMapImpl::remove_short_polylines_out_from_hubs() {
	const ulong polylines_amount = polylines.size();
	for(ulong i=0; i < polylines_amount; i++)
		if(is_valid_polyline(i)) {
			const PolyLineGraph2D::polyline &p = polylines[i];
			if(p.length < SHORT_POLYLINES_MAXLENGTH && (is_extreme(p.start) || is_extreme(p.end)))
				remove_polyline(i);
		}
}

#define MINSPLITLOOP_LENGTH 10

void PolyLineGraph2DHMapImpl::split_loop(ulong pl_id) {
	const PolyLineGraph2D::polyline &p = polylines[pl_id];
	bool tmp;
	const PolyLineGraph2D::polyline::pl_point midpoint = p.next_pl_point_by_length(p.get_extreme_plp(p.start),p.end,p.length/2,tmp);
	if(!tmp)
		split_polyline(PolyLineGraph2D::plg_point(pl_id,midpoint));
}

void PolyLineGraph2DHMapImpl::split_loops() {
	const ulong polylines_amount = polylines.size();
	for(ulong i=0; i < polylines_amount; i++)
		if(is_valid_polyline(i)) {
			const PolyLineGraph2D::polyline &p = polylines[i];
			if(p.length >= MINSPLITLOOP_LENGTH && p.is_loop())
				split_loop(i);
		}
}

void PolyLineGraph2DHMapImpl::optimize() {

	remove_invalid_polylines();
	remove_degenerate_loops();
	remove_2connection_nodes();
	PolyLineGraph2D::optimize();
	connect_close_extremes();
	PolyLineGraph2D::optimize();
	split_loops();

	filter_components_by_polylinesmoothlength(); // 95%
}

void serialize_plg(const PolyLineGraph2DHMapImpl &plg, const string plg_path) {
	std::ofstream ofs(plg_path);
    boost::archive::text_oarchive oa(ofs);
    oa << plg;
}

PolyLineGraph2DHMapImpl deserialize_plg(const string plg_path) {
	PolyLineGraph2DHMapImpl plg;
    std::ifstream ifs(plg_path);
    boost::archive::text_iarchive ia(ifs);
    ia >> plg;
    return plg;
}

void PolyLineGraph2DHMapImpl::prolong_polyline_extreme_and_intersect(const ulong polyline_id, ulong extreme_id, float max_dist)
{
	polyline &pl = polylines[polyline_id];
	pair<glm::vec2,float> p = pl.get_extreme_direction_length_given_length(extreme_id,PROLONG_EXTREME_MIN_SEGMENT_LENGTH);

	if(p.second >= PROLONG_EXTREME_MIN_SEGMENT_LENGTH) {
		plg_point plgp;
		float dist;
		bool found;
		ray2d r(nodes_coords[extreme_id],p.first);
		intersect_ray_first_polyline_within_dist(r, polyline_id, max_dist, plgp, dist, found);
		if(found) {
			cout << "Creating new connection...\n";

			// Split intersected polyline
			ulong split_id = split_polyline(plgp);

			// Connect extreme to split_id
			add_direct_connection(extreme_id,split_id);
		}
	}
}

void PolyLineGraph2DHMapImpl::prolong_polyline_extremes_and_intersect(const ulong polyline_id, float max_dist)
{
	if(!is_valid_polyline(polyline_id))
		return;

	if(is_extreme(polylines[polyline_id].start))
		prolong_polyline_extreme_and_intersect(polyline_id,polylines[polyline_id].start,max_dist);
	if(is_extreme(polylines[polyline_id].end))
		prolong_polyline_extreme_and_intersect(polyline_id,polylines[polyline_id].end,max_dist);
}

void PolyLineGraph2DHMapImpl::prolong_extremes_and_intersect(float max_dist)
{
	cout << "Prolonging polyline extremes...\n";
	for(ulong i=0; i < polylines.size();i++)
		prolong_polyline_extremes_and_intersect(i,max_dist);
}

void PolyLineGraph2DHMapImpl::connect_close_extremes_following_direction() {
	vector<pair<ulong,ulong>> extreme_pairs = find_closest_pairs_with_max_dist_following_direction(get_extreme_nodes_ids_and_coords_and_direction(),DIRECT_CONNECTION_EXTREMES_FOLLOWING_DIRECTION_MAXDIST,DIRECT_CONNECTION_EXTREMES_FOLLOWING_DIRECTION_MINCOS);

	pair<vector<ulong>,vector<set<ulong>>> components = compute_components();
	vector<ulong> node_components_ids = components.first;
	vector<set<ulong>> components_nodes = components.second;

	ulong component_to_change;
	ulong new_id;

	for(const auto &pp : extreme_pairs)
		if(node_components_ids[std::get<0>(pp)] != node_components_ids[std::get<1>(pp)]) {
			if(intersect_polylines(glm::vec4(nodes_coords[std::get<0>(pp)],nodes_coords[std::get<1>(pp)])).size() == 0) {
				add_direct_connection(std::get<0>(pp),std::get<1>(pp));
				if(components_nodes[node_components_ids[std::get<0>(pp)]].size() < components_nodes[node_components_ids[std::get<1>(pp)]].size()) {
					new_id = node_components_ids[std::get<1>(pp)];
					component_to_change = node_components_ids[std::get<0>(pp)];
				} else {
					new_id = node_components_ids[std::get<0>(pp)];
					component_to_change = node_components_ids[std::get<1>(pp)];
				}
				set<ulong> &to_change = components_nodes[component_to_change];
				for(std::set<ulong>::iterator it = to_change.begin(); it != to_change.end(); it++) {
					node_components_ids[*it] = new_id;
				}
			}
		}
}
