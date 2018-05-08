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

#include "glm.hpp"

#include <algorithm>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include "polyline_graph_3d_hmap_impl.hpp"

PolyLineGraph3DHMapImpl::PolyLineGraph3DHMapImpl() {
}

PolyLineGraph3DHMapImpl::PolyLineGraph3DHMapImpl(const vector<polyline> &polylines, const vector<vector<ulong>> &connections, const vector<bool> &visited_nodes, const ulong &real_nodes_amount, const ulong &nodes_amount, const ulong &next_node_id, const vector<glm::vec3> &nodes_coords, const pointmap3dtype &point_map) :
	PolyLineGraph3D(polylines, connections, visited_nodes, real_nodes_amount, nodes_amount, next_node_id, nodes_coords), point_map(point_map) {}

PolyLineGraph3DHMapImpl::~PolyLineGraph3DHMapImpl() {}

ulong PolyLineGraph3DHMapImpl::get_node_id(const glm::vec3 &p_coords) {
	pointmap3dtype::const_iterator it = point_map.find(p_coords);
	ulong node_id;

	if(it != point_map.end() && !is_valid_node(it->second)) {
		invalidate_node(it->second);
		it = point_map.end();
	}

	if(it == point_map.end()) {
		// Create new node, since the coordinates are new
		node_id = get_next_node_id();
		point_map[p_coords] = node_id;
		connections.push_back(vector<ulong>());
		observations.push_back(pair<vector<glm::vec2>,vector<int>>());
		PolyLineGraph3D::add_node_coords(p_coords);
		increment_real_nodes_amount();
	} else
		node_id = it->second;

	return node_id;
}

void PolyLineGraph3DHMapImpl::remap_node(const ulong node_id, const glm::vec3 &new_coords) {
	glm::vec3 old_coords = nodes_coords[node_id];
	remap_node(old_coords,new_coords);
}

void PolyLineGraph3DHMapImpl::remap_node(const glm::vec3 &old_coords, const glm::vec3 &new_coords) {
	pointmap3dtype::const_iterator it = point_map.find(old_coords);
	if(it == point_map.end())
		std::invalid_argument("remap_node - re-mapping inexistent node");
	ulong node_id = it->second;
	point_map.erase(it);
	nodes_coords[node_id]=new_coords;
	point_map[new_coords]=node_id;
}

void PolyLineGraph3DHMapImpl::internal_add_polyline(const polyline &pl) {
	if(!is_duplicate(pl)) {
		ulong pl_id = polylines.size();
		polylines.push_back(pl);
		connections[pl.start].push_back(pl_id);
		if(pl.start != pl.end)
			connections[pl.end].push_back(pl_id);
	}
}

/**
 * Note: check if identical polyline is already in
 */
void PolyLineGraph3DHMapImpl::add_polyline(const vector<glm::vec3> &polyline_coords) {
	//const vector<glm::vec3> polyline_coords2 = filter_polyline(polyline_coords);
	const vector<glm::vec3> &polyline_coords2 = polyline_coords;

	ulong start = get_node_id(polyline_coords2[0]);
	ulong end = get_node_id(polyline_coords2[polyline_coords2.size()-1]);
	polyline pl(start,end,polyline_coords2);
	internal_add_polyline(pl);
}

bool PolyLineGraph3DHMapImpl::is_duplicate(const polyline &pl) {
	const vector<ulong> &s_connections = connections[pl.start];
	const vector<ulong> &e_connections = connections[pl.end];

	const vector<ulong> &smallest_connections = s_connections.size() < e_connections.size() ? s_connections : e_connections;

	for(const auto pl_id : smallest_connections)
		if(polylines[pl_id] == pl)
			return true;

	return false;
}

void PolyLineGraph3DHMapImpl::add_direct_connection(const ulong start, const ulong end) {
	vector<glm::vec3> connection;
	connection.push_back(nodes_coords[start]);
	connection.push_back(nodes_coords[end]);
	polyline pl(start,end,connection);
	internal_add_polyline(pl);
}

void PolyLineGraph3DHMapImpl::add_direct_connection(const glm::vec3 &start_coords, const glm::vec3 &end_coords) {
	ulong start = get_node_id(start_coords);
	ulong end = get_node_id(end_coords);
	add_direct_connection(start,end);
}

void PolyLineGraph3DHMapImpl::add_direct_connection(const glm::vec3 &start_coords, const glm::vec3 &end_coords, pair<ulong, ulong> &out_ids) {
	ulong start = get_node_id(start_coords);
	ulong end = get_node_id(end_coords);
	out_ids.first = start;
	out_ids.second = end;
	add_direct_connection(start,end);
}

void PolyLineGraph3DHMapImpl::remove_polylines_with_longsegments(const float toplength_ratio) {
	vector<float> maxlengths;
	for(ulong pl_id=0;pl_id<polylines.size();pl_id++)
		if(is_valid_polyline(pl_id))
			maxlengths.push_back(polylines[pl_id].get_maxlength());

	const ulong maxlengths_filter_index = maxlengths.size() * toplength_ratio;
	std::nth_element (maxlengths.begin(), maxlengths.begin()+maxlengths_filter_index, maxlengths.end());
	float length_filter = maxlengths[maxlengths_filter_index];

	for(ulong pl_id=0;pl_id<polylines.size();pl_id++)
		if(is_valid_polyline(pl_id) && polylines[pl_id].get_maxlength() >= length_filter)
			remove_polyline(pl_id);
}

void PolyLineGraph3DHMapImpl::filter_nodes(const set3dpoints &inliers) {
	ulong initnodes=nodes_coords.size();
	for(ulong node_id=0;node_id<initnodes;node_id++)
		if(is_valid_node(node_id) && (inliers.find(nodes_coords[node_id])==inliers.end()))
			// node is not an inlier, must be removed
			invalidate_node(node_id);
}

void PolyLineGraph3DHMapImpl::invalidate_node(ulong node_id) {
	point_map.erase(nodes_coords[node_id]);
	PolyLineGraph3D::invalidate_node(node_id);
}

void PolyLineGraph3DHMapImpl::remove_invalid_polylines() {
	const ulong polylines_amount = polylines.size();
	for(ulong i=0; i < polylines_amount; i++)
		if(!is_valid_polyline(i))
			remove_polyline(i);

}

#define SHORT_POLYLINES_MAXLENGTH 2.0

void serialize_plg(const PolyLineGraph3DHMapImpl &plg, const string plg_path) {
	std::ofstream ofs(plg_path);
    boost::archive::text_oarchive oa(ofs);
    oa << plg;
}

PolyLineGraph3DHMapImpl deserialize_3dplg(const string plg_path) {
	PolyLineGraph3DHMapImpl plg;
    std::ifstream ifs(plg_path);
    boost::archive::text_iarchive ia(ifs);
    ia >> plg;
    return plg;
}
