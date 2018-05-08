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

#include <algorithm>
#include <iostream>
#include <iterator>
#include <set>
#include <stack>
#include <stdexcept>

#include "polyline_graph_3d.hpp"

#include "geometric_utilities.hpp"
#include "glm.hpp"


PolyLineGraph3D::PolyLineGraph3D() : nodes_amount(0), next_node_id(0), real_nodes_amount(0) {}

void PolyLineGraph3D::increment_real_nodes_amount() {
	real_nodes_amount++;
}

void PolyLineGraph3D::decrement_real_nodes_amount() {
	real_nodes_amount--;
}

ulong PolyLineGraph3D::get_real_nodes_amount() {
	return real_nodes_amount;
}

void PolyLineGraph3D::polyline::update_length() {
	length = 0.0;
	for(ulong i=1; i < polyline_coords.size();i++)
		length += compute_3d_distance(polyline_coords[i],polyline_coords[i-1]);
}

bool PolyLineGraph3D::polyline::is_start(const ulong node_id) {
	return node_id == start;
}

bool PolyLineGraph3D::polyline::is_end(const ulong node_id) {
	return node_id == end;
}

PolyLineGraph3D::polyline::pl_point::pl_point(const ulong segment_index, const glm::vec3 &coords) : segment_index(segment_index), coords(coords) {}

PolyLineGraph3D::polyline::polyline() {}

PolyLineGraph3D::polyline::polyline(const ulong s, const ulong e, const vector<glm::vec3> &pcs) : start(s), end(e), polyline_coords(pcs) {
	update_length();
}

PolyLineGraph3D::polyline::polyline(const ulong &s, const ulong &e, const vector<glm::vec3> &pcs, float length) : start(s), end(e), polyline_coords(pcs), length(length) { }

bool PolyLineGraph3D::polyline::connects(const ulong a) {
	return a == start || a == end;
}

bool PolyLineGraph3D::polyline::connects(const ulong a, const ulong b) {
	return (a == start && b == end) || (b == start && a == end);
}

float PolyLineGraph3D::polyline::get_maxlength() {
	float max_length=0;
	float curlen;
	for(int i=0; i < polyline_coords.size()-1;i++) {
		curlen = compute_3d_distance(polyline_coords[i],polyline_coords[i+1]);
		max_length = max_length < curlen ? curlen : max_length;
	}
	return max_length;
}

void PolyLineGraph3D::polyline::fragment(const float maxlen) {
	vector<glm::vec3> new_polyline_coords;
	glm::vec3 cur_coords = polyline_coords[0];

	float curlen,nextlen,ratio;
	new_polyline_coords.push_back(cur_coords);
	int cursegm=0;
	while(true) {
		curlen=0;
		while(cursegm<polyline_coords.size()-1 && curlen+compute_3d_distance(cur_coords,polyline_coords[cursegm+1])<maxlen) {
			curlen+=compute_3d_distance(polyline_coords[cursegm],polyline_coords[cursegm+1]);
			cursegm++;
			cur_coords=polyline_coords[cursegm];
		}
		if(cursegm<polyline_coords.size()-1) {
			nextlen=compute_3d_distance(cur_coords,polyline_coords[cursegm+1]);
			ratio=(maxlen-curlen)/(nextlen-curlen);
			cur_coords=cur_coords*(1-ratio)+polyline_coords[cursegm+1]*ratio;
			new_polyline_coords.push_back(cur_coords);
		} else
			break;
	}
	new_polyline_coords.push_back(polyline_coords[polyline_coords.size()-1]);
	polyline_coords=new_polyline_coords;
}

PolyLineGraph3D::plg_point::plg_point() {}

PolyLineGraph3D::plg_point::plg_point(const ulong input_polyline_id, const ulong input_segment_index, const glm::vec3 input_coords) :
		polyline_id(input_polyline_id), segment_index(input_segment_index), coords(input_coords) {}

PolyLineGraph3D::plg_point::plg_point(const pair<ulong,ulong> &input_polyline_segment_index, const glm::vec3 input_coords) :
				polyline_id(input_polyline_segment_index.first), segment_index(input_polyline_segment_index.second), coords(input_coords) {}

vector<vec6> PolyLineGraph3D::polyline::get_segments_list() {
	vector<vec6> res;

	for(ulong i = 1; i < polyline_coords.size(); i++) {
		res.push_back(vec6(polyline_coords[i-1],polyline_coords[i]));
	}

	return res;
}

ulong PolyLineGraph3D::polyline::get_other_end(ulong extreme) {
	return extreme != start ? start : end;
}

bool linearizable_polyline(const vector<glm::vec3> &polyline_coords, const ulong start, const ulong end, const float max_linearizability_distsq) {
	vec6 line = compute_3dline(polyline_coords[start],polyline_coords[end]);

	for(ulong i = start+1; i < end; i++)
		if(distance_point_line_sq(polyline_coords[i], line) > max_linearizability_distsq)
			// current ppline is not compatible with computed 2D line
			return false;

	return true;
}

ulong find_max_se(const vector<glm::vec3> &polyline_coords, const ulong start, const ulong max_se, const float max_linearizability_distsq) {
	if(max_se <= start)
		return start;
	for(ulong cur_se = max_se; cur_se > start+1; cur_se--)
		if(linearizable_polyline(polyline_coords,start,cur_se,max_linearizability_distsq))
			return cur_se;
	return start+1;
}

ulong find_min_eb(const vector<glm::vec3> &polyline_coords, const ulong end, const ulong min_eb, const float max_linearizability_distsq) {
	if(min_eb >= end)
		return end;
	for(ulong cur_eb = min_eb; cur_eb < end - 1; cur_eb++)
		if(linearizable_polyline(polyline_coords,cur_eb,end,max_linearizability_distsq))
			return cur_eb;
	return end-1;
}

/*
 * if start >= end, returns (start,start)
 */
pair<ulong,ulong> find_compatible_se_eb(const vector<glm::vec3> &polyline_coords, const ulong start, const ulong end, const float max_linearizability_distsq) {
	ulong se, eb;
	ulong max_se, min_eb;

	if(start >= end)
		return make_pair(start,start);

	max_se = end;
	min_eb = start;

	do {
		se = find_max_se(polyline_coords, start, max_se,max_linearizability_distsq);
		if(se == end) {
			// The interval start - end is fully linearizable
			eb = 0;
			break;
		}

		eb = find_min_eb(polyline_coords,end,min_eb,max_linearizability_distsq);

		max_se--; // For next iteration, if necessary
		min_eb++;  // For next iteration, if necessary
	} while(eb < se);
	return make_pair(se,eb);
}

vector<glm::vec3> simplify_polyline(const vector<glm::vec3> &polyline_coords, const float max_linearizability_dist) {
	ulong start,end;
	pair<ulong,ulong> se_eb;
	ulong se, eb;
	ulong max_se, min_eb;
	vector<glm::vec3> simplified_polyline_coords, simplified_polyline_coords_end;

	const float max_linearizability_distsq = max_linearizability_dist * max_linearizability_dist;

	start = 0;
	end = polyline_coords.size() - 1;

	simplified_polyline_coords.push_back(polyline_coords[start]);
	simplified_polyline_coords_end.push_back(polyline_coords[end]);

	while(end > start+1) {
		//cout << "Symplifying interval : " << polyline_coords[start] << " - " << polyline_coords[end] << endl;

		se_eb = find_compatible_se_eb(polyline_coords,start,end, max_linearizability_distsq);
		se = se_eb.first;
		eb = se_eb.second;

		if(se == end) {
			// whole interval is linearizable -> DONE
			break;

		} else {
			if( se == eb) {
				// interval is split in: ****** (se = eb) ******
				simplified_polyline_coords.push_back(polyline_coords[se]);
			} else {
				// interval is split in: ****** se ****** eb ****** (with a variable number of elements in the middle, potentially zero)
				simplified_polyline_coords.push_back(polyline_coords[se]);
				simplified_polyline_coords_end.push_back(polyline_coords[eb]);
			}
		}

		start = se;
		end = eb;
	}

	for(vector<glm::vec3>::reverse_iterator it = simplified_polyline_coords_end.rbegin(); it != simplified_polyline_coords_end.rend(); it++)
		simplified_polyline_coords.push_back(*it);

	return simplified_polyline_coords;
}

void PolyLineGraph3D::polyline::simplify() {
	//cout << "Symplifying polyline : " << start << " - " << end << endl;
	simplify(MAXIMUM_LINEARIZABILITY_DISTANCE);
}

void PolyLineGraph3D::polyline::simplify(const float max_linearizability_dist) {
	vector<glm::vec3> spl = simplify_polyline(polyline_coords,max_linearizability_dist);
	polyline_coords.clear();
	polyline_coords = spl;
	update_length();
}

bool PolyLineGraph3D::polyline::operator==(const polyline& p) const
{
    return (start == p.start && end == p.end && vec_glm_vec3_equal(polyline_coords,p.polyline_coords)) || (start == p.end && end == p.start && vec_glm_vec3_equal_inv(polyline_coords,p.polyline_coords));
}

bool PolyLineGraph3D::polyline::is_loop() {
	return start == end;
}

void PolyLineGraph3D::polyline::invalidate() {
	clear_coords();
	length = INVALID_POLYLINE_LENGTH;
}

void PolyLineGraph3D::remove_connection(const ulong node_id, const ulong polyline_id) {
	connections[node_id].erase( std::remove( connections[node_id].begin(), connections[node_id].end(), polyline_id ), connections[node_id].end() );

	if(connections[node_id].size() == 0)
		invalidate_node(node_id);
}

void PolyLineGraph3D::remove_polyline(ulong polyline_id) {
	PolyLineGraph3D::polyline &p = polylines[polyline_id];
	remove_connection(p.start,polyline_id);
	remove_connection(p.end,polyline_id);
	p.invalidate();
}

void PolyLineGraph3D::polyline::clear_coords() {
	polyline_coords.clear();
	polyline_coords.resize(0);
	update_length();
}

PolyLineGraph3D::polyline PolyLineGraph3D::polyline::merge_polylines(const polyline& p1,const polyline& p2) {
	ulong common_extreme;
	ulong p3s,p3e;
	vector<glm::vec3> p3coords;

	if(p1.start == p2.start) {
		//  <---- P1 P2 ---->
		p3s = p1.end;
		p3e = p2.end;
		for(std::vector<glm::vec3>::const_reverse_iterator rit = p1.polyline_coords.rbegin(); rit!= p1.polyline_coords.rend(); ++rit)
			p3coords.push_back(*rit);
		std::vector<glm::vec3>::const_iterator it = p2.polyline_coords.begin();
		it++;
		for (; it != p2.polyline_coords.end(); ++it)
			p3coords.push_back(*it);
		return polyline(p3s,p3e,p3coords);
	} else if(p1.start == p2.end) {
		//  P2 ----> P1 ---->

		p3s = p2.start;
		p3e = p1.end;
		for(std::vector<glm::vec3>::const_iterator it = p2.polyline_coords.begin(); it!= p2.polyline_coords.end(); ++it)
			p3coords.push_back(*it);
		std::vector<glm::vec3>::const_iterator it = p1.polyline_coords.begin();
		it++;
		for (; it != p1.polyline_coords.end(); ++it)
			p3coords.push_back(*it);
		return polyline(p3s,p3e,p3coords);

	} else if(p1.end == p2.start) {
		//  P1 ----> P2 ---->

		p3s = p1.start;
		p3e = p2.end;
		for(std::vector<glm::vec3>::const_iterator it = p1.polyline_coords.begin(); it!= p1.polyline_coords.end(); ++it)
			p3coords.push_back(*it);
		std::vector<glm::vec3>::const_iterator it = p2.polyline_coords.begin();
		it++;
		for (; it != p2.polyline_coords.end(); ++it)
			p3coords.push_back(*it);
		return polyline(p3s,p3e,p3coords);

	} else if(p1.end == p2.end) {
		//  P1 ----> <---- P2
		p3s = p1.start;
		p3e = p2.start;
		for(std::vector<glm::vec3>::const_iterator it = p1.polyline_coords.begin(); it!= p1.polyline_coords.end(); ++it)
			p3coords.push_back(*it);
		std::vector<glm::vec3>::const_reverse_iterator rit = p2.polyline_coords.rbegin();
		rit++;
		for (; rit != p2.polyline_coords.rend(); ++rit)
			p3coords.push_back(*rit);
		return polyline(p3s,p3e,p3coords);
	} else
		throw std::invalid_argument( "cannot merge disconnected polylines" );
}

void PolyLineGraph3D::simplify()
{
	simplify(MAXIMUM_LINEARIZABILITY_DISTANCE);
}

void PolyLineGraph3D::simplify(const float max_linearizability_dist)
{
	for(ulong pl_id=0; pl_id < polylines.size(); pl_id++)
		if(is_valid_polyline(pl_id))
			polylines[pl_id].simplify(max_linearizability_dist);
}

void PolyLineGraph3D::invalidate_node(ulong node_id) {
	nodes_coords[node_id] = glm::vec3(INVALID_POINT_COORDS,INVALID_POINT_COORDS,INVALID_POINT_COORDS);
	connections[node_id].clear();
	connections[node_id].resize(0);
}

bool PolyLineGraph3D::is_valid_node(const ulong node_id) {
	return nodes_coords[node_id][0] != INVALID_POINT_COORDS && nodes_coords[node_id][1] != INVALID_POINT_COORDS;
}

bool PolyLineGraph3D::is_valid_polyline(const ulong polyline_id) {
	return is_valid_node(polylines[polyline_id].start)
			&& is_valid_node(polylines[polyline_id].end)
			&& polylines[polyline_id].polyline_coords.size() > 1
			&& get_node_coords(polylines[polyline_id].start) == polylines[polyline_id].polyline_coords[0]
			&& get_node_coords(polylines[polyline_id].end) == polylines[polyline_id].polyline_coords[polylines[polyline_id].polyline_coords.size()-1];
}

vector<vec6> PolyLineGraph3D::get_segments_list() {
	vector<vec6> res;

	for(ulong i=0; i < polylines.size(); i++)
		if(is_valid_polyline(i)) {
			const polyline &p = polylines[i];
			const vector<vec6> cur_res = p.get_segments_list();
			for(const auto &s : cur_res)
				res.push_back(s);
		}

	return res;
}

vector<vector<vec6>> PolyLineGraph3D::get_segments_list_by_polyline() {
	vector<vector<vec6>> res;

	for(ulong i=0; i < polylines.size(); i++)
		if(is_valid_polyline(i)) {
			const polyline &p = polylines[i];
			res.push_back(p.get_segments_list());
		}

	return res;
}

vector<vec6> PolyLineGraph3D::get_segments_list(set3dpoints include_only) {
	vector<vec6> res;

	for(ulong i=0; i < polylines.size(); i++)
		if(is_valid_polyline(i)) {
			const polyline &p = polylines[i];
			const vector<vec6> cur_res = p.get_segments_list();
			for(const auto &s : cur_res)
				if((include_only.find(s.start()) != include_only.end()) && (include_only.find(s.end()) != include_only.end()))
					res.push_back(s);
		}

	return res;
}

glm::vec3 PolyLineGraph3D::get_node_coords(ulong node_id) {
	return nodes_coords[node_id];
}

vector<glm::vec3> PolyLineGraph3D::get_nodes_coords() {
	return nodes_coords;
}

ulong PolyLineGraph3D::get_nodes_amount() {
	return nodes_amount;
}

void PolyLineGraph3D::fragment(const float maxlen) {
	for(ulong pl_id=0; pl_id < polylines.size(); pl_id++)
		if(is_valid_polyline(pl_id))
			polylines[pl_id].fragment(maxlen);
}

ulong PolyLineGraph3D::get_polylines_amount() {
	return polylines.size();
}

void PolyLineGraph3D::add_node_coords(const glm::vec3 &p_coords)
{
	nodes_coords.push_back(p_coords);
}

vector<vector<ulong>> PolyLineGraph3D::get_connections() {
	return connections;
}

vector<PolyLineGraph3D::polyline> PolyLineGraph3D::get_polylines() {
	return polylines;
}

ulong PolyLineGraph3D::get_next_node_id() {
	ulong cur = next_node_id;
	next_node_id++;
	nodes_amount++;
	visited_nodes.push_back(false);
	return cur;
}

bool PolyLineGraph3D::has_loop(ulong node_id)
{
	for(const auto connection : connections[node_id])
		if(polylines[connection].is_loop())
			return true;
	return false;
}

void PolyLineGraph3D::set_observations(ulong node_id, vector<glm::vec2> projections, vector<int> cam_ids) {
	observations[node_id].first = projections;
	observations[node_id].second = cam_ids;
}

/* Hub, either:
 * - >= 3 connected polylines
 * - 2 connected polylines, provided 1 polyline is a loop
 */
bool PolyLineGraph3D::is_hub(ulong node_id)
{
	const int connections_amount = connections[node_id].size();
	return (connections_amount > 2) || (connections_amount == 2 && has_loop(node_id));
}

/* Extreme, either:
 * - 1 connected polyline (not a loop)
 */
bool PolyLineGraph3D::is_extreme(ulong node_id)
{
	return connections[node_id].size() == 1 && !polylines[connections[node_id][0]].is_loop();
}

/* LoopNode, either:
 * - 1 connected polyline (a loop)
 */
bool PolyLineGraph3D::is_loopnode(ulong node_id)
{
	return connections[node_id].size() == 1 && polylines[connections[node_id][0]].is_loop();
}

vector<ulong> PolyLineGraph3D::get_hub_nodes()
{
	vector<ulong> res;
	for(ulong node_id = 0; node_id < get_nodes_amount(); node_id++)
		if(is_valid_node(node_id) && is_hub(node_id))
			res.push_back(node_id);
	return res;
}

vector<ulong> PolyLineGraph3D::get_extreme_nodes()
{
	vector<ulong> res;
	for(ulong node_id = 0; node_id < get_nodes_amount(); node_id++)
		if(is_valid_node(node_id) && is_extreme(node_id))
			res.push_back(node_id);
	return res;
}

vector<ulong> PolyLineGraph3D::get_loopnodes()
{
	vector<ulong> res;
	for(ulong node_id = 0; node_id < get_nodes_amount(); node_id++)
		if(is_valid_node(node_id) && is_loopnode(node_id))
			res.push_back(node_id);
	return res;
}

vector<ulong> PolyLineGraph3D::get_nodes_with_loops()
{
	vector<ulong> res;
	for(ulong node_id = 0; node_id < get_nodes_amount(); node_id++)
		if(is_valid_node(node_id) && has_loop(node_id))
			res.push_back(node_id);
	return res;
}

vector<glm::vec3> PolyLineGraph3D::get_hub_nodes_coords()
{
	vector<glm::vec3> res;
	for(ulong node_id = 0; node_id < get_nodes_amount(); node_id++)
		if(is_valid_node(node_id) && is_hub(node_id))
			res.push_back(nodes_coords[node_id]);
	return res;
}

vector<glm::vec3> PolyLineGraph3D::get_extreme_nodes_coords()
{
	vector<glm::vec3> res;
	for(ulong node_id = 0; node_id < get_nodes_amount(); node_id++)
		if(is_valid_node(node_id) && is_extreme(node_id))
			res.push_back(nodes_coords[node_id]);
	return res;
}

vector<glm::vec3> PolyLineGraph3D::get_loopnodes_coords()
{
	vector<glm::vec3> res;
	for(ulong node_id = 0; node_id < get_nodes_amount(); node_id++)
		if(is_valid_node(node_id) && is_loopnode(node_id))
			res.push_back(nodes_coords[node_id]);
	return res;
}

vector<glm::vec3> PolyLineGraph3D::get_nodes_with_loops_coords()
{
	vector<glm::vec3> res;
	for(ulong node_id = 0; node_id < get_nodes_amount(); node_id++)
		if(is_valid_node(node_id) && has_loop(node_id))
			res.push_back(nodes_coords[node_id]);
	return res;
}

bool PolyLineGraph3D::is_connected_node(const ulong start_node,const ulong end_node) {
	vector<bool> visited_vec;
	visited_vec.push_back(start_node);
	visited_nodes[start_node]=true;

	bool found = false;

	stack<ulong> cur_to_visit;
	cur_to_visit.push(start_node);
	ulong cur_dist = 0;

	while(!found && !cur_to_visit.empty()) {
		const ulong cur_node = cur_to_visit.top();
		cur_to_visit.pop();

		for(const auto cur_polyline_id : connections[cur_node]) {
			const ulong connected_node = polylines[cur_polyline_id].get_other_end(cur_node);

			if(connected_node == end_node) {
				found = true;
				break;
			}

			if(!visited_nodes[connected_node]) {
				cur_to_visit.push(connected_node);
				visited_nodes[connected_node] = true;
				visited_vec.push_back(connected_node);
			}
		}
	}


	for(auto v : visited_vec)
		visited_nodes[v] = false;

	return found;
}

PolyLineGraph3D::PolyLineGraph3D(const vector<polyline> &polylines, const vector<vector<ulong>> &connections, const vector<bool> &visited_nodes, const ulong &real_nodes_amount, const ulong &nodes_amount, const ulong &next_node_id, const vector<glm::vec3> &nodes_coords) :
		polylines(polylines), connections(connections), visited_nodes(visited_nodes), real_nodes_amount(real_nodes_amount), nodes_amount(nodes_amount), next_node_id(next_node_id), nodes_coords(nodes_coords) {}

bool PolyLineGraph3D::is_connected_node(ulong start_node,ulong end_node, ulong max_jumps) {
	vector<bool> visited_vec;
	visited_vec.push_back(start_node);
	visited_nodes[start_node]=true;

	bool found = false;

	stack<ulong> cur_to_visit;
	cur_to_visit.push(start_node);
	ulong cur_dist = 0;

	while(cur_dist <= max_jumps && !cur_to_visit.empty()) {
		stack<ulong> next_to_visit;

		while(!found && !cur_to_visit.empty()) {
			const ulong cur_node = cur_to_visit.top();
			cur_to_visit.pop();

			for(const auto cur_polyline_id : connections[cur_node]) {
				const ulong connected_node = polylines[cur_polyline_id].get_other_end(cur_node);

				if(connected_node == end_node) {
					found = true;
					break;
				}

				if(!visited_nodes[connected_node]) {
					next_to_visit.push(connected_node);
					visited_nodes[connected_node] = true;
					visited_vec.push_back(connected_node);
				}
			}
		}

		cur_to_visit = next_to_visit;
		cur_dist++;
	}

	for(auto v : visited_vec)
		visited_nodes[v] = false;

	return found;
}

bool can_follow_polyline_without_going_outside_radius(const vector<glm::vec3> &polyline_coords,const glm::vec3 &a,const float detection_radius_sq) {
	for(const auto &p : polyline_coords)
		if(squared_3d_distance(p,a) > detection_radius_sq)
			return false;

	return true;
}

bool PolyLineGraph3D::is_connected_node_inside_radius(const ulong start_node,const ulong end_node, const float detection_radius) {
	vector<bool> visited_vec;
	visited_vec.push_back(start_node);
	visited_nodes[start_node]=true;

	const glm::vec3 &start_coords = nodes_coords[start_node];

	const float detection_radius_sq = detection_radius*detection_radius;

	if(squared_3d_distance(nodes_coords[start_node],nodes_coords[end_node]) > detection_radius_sq)
		return false;

	bool found = false;

	stack<ulong> cur_to_visit;
	cur_to_visit.push(start_node);
	ulong cur_dist = 0;

	while(!found && !cur_to_visit.empty()) {
		const ulong cur_node = cur_to_visit.top();
		cur_to_visit.pop();

		for(const auto cur_polyline_id : connections[cur_node]) {
			const ulong connected_node = polylines[cur_polyline_id].get_other_end(cur_node);

			if(connected_node == end_node) {
				found = true;
				break;
			}

			if(!visited_nodes[connected_node] && can_follow_polyline_without_going_outside_radius(polylines[cur_polyline_id].polyline_coords,start_coords,detection_radius_sq)) {
				cur_to_visit.push(connected_node);
				visited_nodes[connected_node] = true;
				visited_vec.push_back(connected_node);
			}
		}
	}


	for(auto v : visited_vec)
		visited_nodes[v] = false;

	return found;
}

bool PolyLineGraph3D::is_connected_polyline(ulong polyline_id_start,ulong polyline_id_end) {
	const polyline &p_start = polylines[polyline_id_start];

	vector<bool> visited_vec;
	visited_vec.push_back(p_start.start);
	visited_vec.push_back(p_start.end);

	visited_nodes[p_start.start]=true;
	visited_nodes[p_start.end]=true;

	bool found = false;

	stack<ulong> cur_to_visit;
	cur_to_visit.push(p_start.start);
	cur_to_visit.push(p_start.end);
	ulong cur_dist = 0;

	while(!found && !cur_to_visit.empty()) {
		const ulong cur_node = cur_to_visit.top();
		cur_to_visit.pop();

		for(const auto cur_polyline_id : connections[cur_node]) {
			if(cur_polyline_id == polyline_id_end) {
				found = true;
				break;
			}

			const ulong connected_node = polylines[cur_polyline_id].get_other_end(cur_node);

			if(!visited_nodes[connected_node]) {
				cur_to_visit.push(connected_node);
				visited_nodes[connected_node] = true;
				visited_vec.push_back(connected_node);
			}
		}
	}


	for(auto v : visited_vec)
		visited_nodes[v] = false;

	return found;
}

bool PolyLineGraph3D::is_connected_polyline(ulong polyline_id_start,ulong polyline_id_end, ulong max_jumps) {
	const polyline &p_start = polylines[polyline_id_start];

	vector<bool> visited_vec;
	visited_vec.push_back(p_start.start);
	visited_vec.push_back(p_start.end);

	visited_nodes[p_start.start]=true;
	visited_nodes[p_start.end]=true;

	bool found = false;

	stack<ulong> cur_to_visit;
	cur_to_visit.push(p_start.start);
	cur_to_visit.push(p_start.end);
	ulong cur_dist = 0;

	while(cur_dist <= max_jumps && !cur_to_visit.empty()) {
		stack<ulong> next_to_visit;

		while(!found && !cur_to_visit.empty()) {
			const ulong cur_node = cur_to_visit.top();
			cur_to_visit.pop();

			for(const auto cur_polyline_id : connections[cur_node]) {
				if(cur_polyline_id == polyline_id_end) {
					found = true;
					break;
				}

				const ulong connected_node = polylines[cur_polyline_id].get_other_end(cur_node);

				if(!visited_nodes[connected_node]) {
					next_to_visit.push(connected_node);
					visited_nodes[connected_node] = true;
					visited_vec.push_back(connected_node);
				}
			}
		}

		for(auto v : visited_vec)
			visited_nodes[v] = false;

		cur_to_visit = next_to_visit;
		cur_dist++;
	}

	return found;
}

// Returns the minimum between distsq(p,x) and distsq(p,y)
inline float min_dist_sq(const glm::vec3 &p,const glm::vec3 &x,const glm::vec3 &y) {
	return min(squared_3d_distance(p,x),squared_3d_distance(p,y));
}

bool can_follow_polyline_without_going_outside_radius(const vector<glm::vec3> &polyline_coords,const glm::vec3 &s_a,const glm::vec3 &s_b,const float detection_radius_sq) {
	for(const auto &p : polyline_coords)
		if(min_dist_sq(p,s_a,s_b) > detection_radius_sq)
			return false;

	return true;
}

bool PolyLineGraph3D::is_connected_polyline_inside_radius(ulong polyline_id_start, ulong polyline_id_end, const float detection_radius) {
	const polyline &p_start = polylines[polyline_id_start];
	const polyline &p_end = polylines[polyline_id_end];

	const float detection_radius_sq = detection_radius*detection_radius;

	const glm::vec3 &s_a = nodes_coords[p_start.start];
	const glm::vec3 &s_b = nodes_coords[p_start.end];
	const glm::vec3 &e_a = nodes_coords[p_end.start];
	const glm::vec3 &e_b = nodes_coords[p_end.end];

	if(min_dist_sq(s_a,e_a,e_b) > detection_radius_sq && min_dist_sq(s_b,e_a,e_b) > detection_radius_sq)
		return false;

	vector<bool> visited_vec;
	visited_vec.push_back(p_start.start);
	visited_vec.push_back(p_start.end);

	visited_nodes[p_start.start]=true;
	visited_nodes[p_start.end]=true;

	bool found = false;

	stack<ulong> cur_to_visit;
	cur_to_visit.push(p_start.start);
	cur_to_visit.push(p_start.end);
	ulong cur_dist = 0;

	while(!found && !cur_to_visit.empty()) {
		const ulong cur_node = cur_to_visit.top();
		cur_to_visit.pop();

		for(const auto cur_polyline_id : connections[cur_node]) {
			if(cur_polyline_id == polyline_id_end) {
				found = true;
				break;
			}

			const ulong connected_node = polylines[cur_polyline_id].get_other_end(cur_node);

			if(!visited_nodes[connected_node] && can_follow_polyline_without_going_outside_radius(polylines[cur_polyline_id].polyline_coords,s_a,s_b,detection_radius_sq)) {
				cur_to_visit.push(connected_node);
				visited_nodes[connected_node] = true;
				visited_vec.push_back(connected_node);
			}
		}
	}


	for(auto v : visited_vec)
		visited_nodes[v] = false;

	return found;
}

void PolyLineGraph3D::print_stats() {
	cout << "PLG 3D STATS" << endl << endl;
	cout << "Nodes: " << nodes_coords.size() << endl;
	cout << "Polylines: " << polylines.size() << endl;
}

vector<ulong> PolyLineGraph3D::valid_nodes_index(ulong &valid_nodes_amount) {
	vector<ulong> res(nodes_amount);
	ulong cur_valid_index=0;

	for(ulong node_id=0; node_id < nodes_amount; node_id++) {
		if(is_valid_node(node_id)) {
			res[node_id] = cur_valid_index;
			cur_valid_index++;
		}
	}

	valid_nodes_amount = cur_valid_index;

	return res;
}

ulong PolyLineGraph3D::get_amount_valid_polylines() {
	ulong cur_valid_index=0;

	for(ulong pl_id=0; pl_id < nodes_amount; pl_id++) {
		if(is_valid_polyline(pl_id)) {
			cur_valid_index++;
		}
	}

	return cur_valid_index;
}
