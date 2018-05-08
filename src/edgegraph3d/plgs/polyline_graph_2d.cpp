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
#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <stack>
#include <stdexcept>

#include "polyline_graph_2d.hpp"

#include "geometric_utilities.hpp"
#include "glm.hpp"


PolyLineGraph2D::PolyLineGraph2D() : nodes_amount(0), next_node_id(0), real_nodes_amount(0) {}

void PolyLineGraph2D::increment_real_nodes_amount() {
	real_nodes_amount++;
}

PolyLineGraph2D::plg_node_components::plg_node_components(vector<ulong> &nodes_component_id,vector<set<ulong>> &components_nodes_ids) :
		nodes_component_id(nodes_component_id), components_nodes_ids(components_nodes_ids) {}

PolyLineGraph2D::plg_node_components::~plg_node_components() {}

PolyLineGraph2D::plg_polyline_components::plg_polyline_components(vector<ulong> &polylines_component_id,vector<set<ulong>> &components_polylines_ids) :
		polylines_component_id(polylines_component_id), components_polylines_ids(components_polylines_ids) {}

PolyLineGraph2D::plg_polyline_components::~plg_polyline_components() {}

PolyLineGraph2D::plg_components::plg_components(plg_node_components &nc, plg_polyline_components &pc) :
		nc(nc), pc(pc) {}

PolyLineGraph2D::plg_components::plg_components(vector<ulong> &nodes_component_id,vector<set<ulong>> &components_nodes_ids, vector<ulong> &polylines_component_id,vector<set<ulong>> &components_polylines_ids) :
		nc(plg_node_components(nodes_component_id,components_nodes_ids)), pc(plg_polyline_components(polylines_component_id, components_polylines_ids)) {}

PolyLineGraph2D::plg_components::~plg_components() {}

void PolyLineGraph2D::decrement_real_nodes_amount() {
	real_nodes_amount--;
}

ulong PolyLineGraph2D::get_real_nodes_amount() {
	return real_nodes_amount;
}

void PolyLineGraph2D::polyline::update_length() {
	length = 0.0;
	for(ulong i=1; i < polyline_coords.size();i++)
		length += compute_2d_distance(polyline_coords[i],polyline_coords[i-1]);
}

float PolyLineGraph2D::polyline::compute_max_smooth_length() {
	float maxl = 0.0;
	float cur_length;

	ulong i=1;
	while(i < polyline_coords.size()) {
		cur_length = compute_2d_distance(polyline_coords[i],polyline_coords[i-1]);
		for(i++; i < polyline_coords.size(); i++)
			if(compute_anglecos_vec2_vec2(polyline_coords[i-1],polyline_coords[i],polyline_coords[i-2],polyline_coords[i-1]))
				cur_length += compute_2d_distance(polyline_coords[i],polyline_coords[i-1]);
			else
				break;

		maxl = maxl < cur_length ? cur_length : maxl;
	}

	return maxl;
}

ulong PolyLineGraph2D::polyline::get_extreme_id(const pl_point &plp) {
	if(is_start(plp))
		return start;
	else if(is_end(plp))
		return end;
	else
		std::invalid_argument("get_extreme_id - nor start or end passed");
}

bool PolyLineGraph2D::polyline::is_start(const ulong node_id) {
	return node_id == start;
}

bool PolyLineGraph2D::polyline::is_start(const pl_point &plp) {
	return plp.segment_index == 0 && plp.coords == polyline_coords[0];
}

bool PolyLineGraph2D::polyline::is_end(const ulong node_id) {
	return node_id == end;
}

bool PolyLineGraph2D::polyline::is_end(const pl_point &plp) {
	return plp.segment_index == polyline_coords.size()-2 && plp.coords == polyline_coords[polyline_coords.size()-1];
}

glm::vec2 PolyLineGraph2D::polyline::get_extreme_coordinates(const ulong extreme_id) {
	if(extreme_id == start)
		return polyline_coords[0];
	else if(extreme_id == end)
		return polyline_coords[polyline_coords.size()-1];
	else
		std::invalid_argument("get_extreme_coordinates - nor start or end passed");
}

PolyLineGraph2D::polyline::pl_point PolyLineGraph2D::polyline::get_start_plp() {
	return pl_point(0,polyline_coords[0]);
}
PolyLineGraph2D::polyline::pl_point PolyLineGraph2D::polyline::get_end_plp() {
	return pl_point(polyline_coords.size()-2,polyline_coords[polyline_coords.size()-1]);
}

PolyLineGraph2D::polyline::pl_point PolyLineGraph2D::polyline::get_extreme_plp(const ulong extreme_id) {
	if(extreme_id == start)
		return get_start_plp();
	else if(extreme_id == end)
		return get_end_plp();
	else
		std::invalid_argument("get_extreme_plp - nor start or end passed");
}

PolyLineGraph2D::polyline::pl_point PolyLineGraph2D::polyline::get_extreme_plp(const ulong extreme_id, bool &valid) {
	valid=false;
	if(extreme_id == start) {
		valid=true;
		return get_start_plp();
	} else if(extreme_id == end) {
		valid=true;
		return get_end_plp();
	}
}

bool PolyLineGraph2D::polyline::has_point(const pl_point &plp) {
	return plp.segment_index < polyline_coords.size() - 1 && aligned(polyline_coords[plp.segment_index],polyline_coords[plp.segment_index+1],plp.coords);
}

pair<vector<glm::vec2>,vector<glm::vec2>> PolyLineGraph2D::polyline::split(const pl_point &plp) {
	if(!has_point(plp))
		std::invalid_argument("polyline::split - plp doesn't belong to polyline");

	std::vector<glm::vec2> coords1;
	for(ulong i=0; i <= plp.segment_index; i++)
		coords1.push_back(polyline_coords[i]);
	if(plp.coords != polyline_coords[plp.segment_index])
		coords1.push_back(plp.coords);

	std::vector<glm::vec2> coords2;
	coords2.push_back(plp.coords);
	for(ulong i=plp.segment_index+1; i < polyline_coords.size(); i++)
		coords2.push_back(polyline_coords[i]);

	return make_pair(coords1,coords2);
}

// Input: start or end id
// Output: direction and length of last segment
pair<glm::vec2,float> PolyLineGraph2D::polyline::get_extreme_direction_length(const ulong extreme) {
	if(is_start(extreme))
		return make_pair(polyline_coords[0]-polyline_coords[1],compute_2d_distance(polyline_coords[0],polyline_coords[1]));
	else if(is_end(extreme)) {
		const ulong sz = polyline_coords.size();
		return make_pair(polyline_coords[sz-1]-polyline_coords[sz-2],compute_2d_distance(polyline_coords[sz-1],polyline_coords[sz-2]));
	} else
		std::invalid_argument("get_extreme_direction_length - nor start or end passed");
}

// Input: start or end id
// Output: direction and length of last portion of polyline with given lengthsq
pair<glm::vec2,float> PolyLineGraph2D::polyline::get_extreme_direction_length_given_length(const ulong extreme, const float length) {
	const ulong sz = polyline_coords.size();
	float cur_lengthsq;
	const float init_lengthsq = length * length;
	float residual_lengthsq = init_lengthsq;
	if(is_start(extreme)) {
		const glm::vec2 init_node = polyline_coords[0];
		glm::vec2 final_node;
		glm::vec2 last_segment;
		ulong i;
		for(i=1; i < sz; i++)
		{
			last_segment = polyline_coords[i] - polyline_coords[i-1];
			cur_lengthsq = compute_lengthsq(last_segment);
			if(residual_lengthsq <= cur_lengthsq)
				break;
			residual_lengthsq -= cur_lengthsq;
		}
		if(i>=sz) {
			// exploration interrupted by reaching other polyline extreme
			final_node = polyline_coords[sz-1];
		} else {
			float ratio = residual_lengthsq / cur_lengthsq;
			residual_lengthsq = 0.0;
			final_node = polyline_coords[i-1] + ratio * last_segment;
		}
		const glm::vec2 dir = init_node - final_node;
		const float followed_length = sqrt(init_lengthsq - residual_lengthsq);

		return make_pair(dir,followed_length);
	} else if(is_end(extreme)) {
		const glm::vec2 init_node = polyline_coords[sz-1];
		glm::vec2 final_node;
		glm::vec2 last_segment;
		ulong i;
		for(i=sz-1; i > 0; i--)
		{
			last_segment = polyline_coords[i-1] - polyline_coords[i];
			cur_lengthsq = compute_lengthsq(last_segment);
			if(residual_lengthsq <= cur_lengthsq)
				break;
			residual_lengthsq -= cur_lengthsq;
		}
		if(i==0) {
			// exploration interrupted by reaching other polyline extreme
			final_node = polyline_coords[0];
		} else {
			float ratio = residual_lengthsq / cur_lengthsq;
			residual_lengthsq = 0.0;
			final_node = polyline_coords[i] + ratio * last_segment;
		}
		const glm::vec2 dir = final_node - init_node;
		const float followed_length = sqrt(init_lengthsq - residual_lengthsq);

		return make_pair(dir,followed_length);
	} else
		std::invalid_argument("get_extreme_direction_length - nor start or end passed");
}

PolyLineGraph2D::polyline::pl_point::pl_point() {}

PolyLineGraph2D::polyline::pl_point::pl_point(const ulong segment_index, const glm::vec2 &coords) : segment_index(segment_index), coords(coords) {}

bool PolyLineGraph2D::polyline::pl_point::operator==(const pl_point& a) const {
	return coords == a.coords && segment_index == a.segment_index;
}

PolyLineGraph2D::polyline::pl_interval::pl_interval() {}

PolyLineGraph2D::polyline::pl_interval::pl_interval(const PolyLineGraph2D::polyline::pl_point start, const PolyLineGraph2D::polyline::pl_point end) : start(start), end(end) {}

bool PolyLineGraph2D::polyline::interval_contains_plp(const pl_interval &pli, const pl_point &plp)  {
	if(plp.segment_index < pli.start.segment_index - 1 || plp.segment_index > pli.end.segment_index + 1)
		return false;
	if(plp.segment_index > pli.start.segment_index && plp.segment_index < pli.start.segment_index)
		return true;
	if(plp.segment_index == pli.start.segment_index) {
		const glm::vec2 &prevp = polyline_coords[pli.start.segment_index];
		if(is_ordered_2dlinepoints(prevp,pli.start.coords,plp.coords)) {
			if(pli.end.segment_index == pli.start.segment_index) {
				return is_ordered_2dlinepoints(pli.start.coords,plp.coords,pli.end.coords);
			} else
				return true;
		} else
			return false;
	}
	if(plp.segment_index == pli.end.segment_index) {
		const glm::vec2 &prevp = polyline_coords[pli.end.segment_index];
		return is_ordered_2dlinepoints(prevp,plp.coords,pli.end.coords);
	}
	if(plp.segment_index == pli.end.segment_index + 1 && polyline_coords[plp.segment_index] == plp.coords && pli.end.coords == polyline_coords[plp.segment_index])
		return true; // coincides with end
	if(plp.segment_index == pli.start.segment_index - 1 && polyline_coords[plp.segment_index+1] == plp.coords && pli.start.coords == polyline_coords[plp.segment_index+1])
		return true; // coincides with start
	return false;
}

vector<PolyLineGraph2D::polyline::pl_point> PolyLineGraph2D::polyline::intersect_segment(const glm::vec4 &segment)
{
	vector<PolyLineGraph2D::polyline::pl_point> res;
	bool parallel, overlapped, intersection_found;
	 glm::vec2 intersection;

	for(ulong i=1; i < polyline_coords.size(); i++)
	{
		intersect_segment_segment(glm::vec4(polyline_coords[i],polyline_coords[i-1]),segment, parallel, overlapped, intersection_found, intersection);
		if(intersection_found) {
			res.push_back(PolyLineGraph2D::polyline::pl_point(i-1,intersection));
		}
	}

	return res;
}

vector<PolyLineGraph2D::polyline::pl_point> PolyLineGraph2D::polyline::intersect_line(const glm::vec3 &line)
{
	vector<PolyLineGraph2D::polyline::pl_point> res;
	bool parallel, overlapped, intersection_found;
	 glm::vec2 intersection;

	for(ulong i=1; i < polyline_coords.size(); i++)
	{
		intersect_segment_line(glm::vec4(polyline_coords[i],polyline_coords[i-1]),line, parallel, overlapped, intersection_found, intersection);
		if(intersection_found) {
			res.push_back(PolyLineGraph2D::polyline::pl_point(i-1,intersection));
		}
	}

	return res;
}

vector<PolyLineGraph2D::polyline::pl_point> PolyLineGraph2D::polyline::intersect_ray(const ray2d &r)
{
	vector<PolyLineGraph2D::polyline::pl_point> res;
	bool parallel, overlapped, past_center, intersection_found;
	 glm::vec2 intersection;

	for(ulong i=1; i < polyline_coords.size(); i++)
	{
		intersect_segment_ray(glm::vec4(polyline_coords[i],polyline_coords[i-1]),r, parallel, overlapped, past_center, intersection_found, intersection);
		if(intersection_found) {
			res.push_back(PolyLineGraph2D::polyline::pl_point(i-1,intersection));
		}
	}

	return res;
}

void PolyLineGraph2D::polyline::first_intersect_ray(const ray2d &r, pl_point &intersection, float &distance, bool &found) {
	vector<PolyLineGraph2D::polyline::pl_point> res = intersect_ray(r);

	if(res.size()==0)
	{
		found = false;
		return;
	}

	int min_ind = -1;
	float min_distsq = std::numeric_limits<float>::max();
	float cur_distsq;
	for(int i=0; i < res.size(); i++) {
		cur_distsq = squared_2d_distance(r.start,res[i].coords);
		if(cur_distsq < min_distsq) {
			cur_distsq = min_distsq;
			min_ind = i;
		}
	}

	distance = sqrt(min_distsq);
	found = true;
	intersection = res[min_ind];
}

void PolyLineGraph2D::polyline::first_intersect_ray_approx(const ray2d &r, pl_point &intersection, float &distance, bool &found) {
	first_intersect_ray(r, intersection, distance, found);

	if(!found)
		return;

	float dprev = squared_2d_distance(intersection.coords,polyline_coords[intersection.segment_index]);
	float dnext = squared_2d_distance(intersection.coords,polyline_coords[intersection.segment_index+1]);

	if(dprev <= dnext && dprev <= POLYLINE_RAY_INTERSECT_APPROX_MAX_DISTSQ)
	{
		intersection = pl_point(intersection.segment_index,polyline_coords[intersection.segment_index]);
		distance = compute_2d_distance(intersection.coords,r.start);
		return;
	} else if (dnext <= POLYLINE_RAY_INTERSECT_APPROX_MAX_DISTSQ) {
		intersection = pl_point(intersection.segment_index,polyline_coords[intersection.segment_index+1]);
		distance = compute_2d_distance(intersection.coords,r.start);
	}
}

PolyLineGraph2D::polyline::pl_point PolyLineGraph2D::polyline::next_pl_point_by_distance(const PolyLineGraph2D::polyline::pl_point init_plp, const ulong direction, const float distance, bool &reached_polyline_extreme)
{
	float prevdist, curdist;
	reached_polyline_extreme = false;
	float ratio;
	ulong i;

	if(direction == start) {
		curdist = compute_2d_distance(polyline_coords[init_plp.segment_index],init_plp.coords);
		if(curdist >= distance) {
			ratio = distance / curdist;
			return PolyLineGraph2D::polyline::pl_point(init_plp.segment_index,first_plus_ratio_of_segment(init_plp.coords,polyline_coords[init_plp.segment_index], ratio));
		}
		for(i=init_plp.segment_index; i > 0; i--) {
			// segment polyline_coords[i-1] <--> polyline_coords[i]
			prevdist = curdist; // Now prevdist is distance between init_plp.coords and polyline_coords[i]
			curdist = compute_2d_distance(polyline_coords[i-1],init_plp.coords); // Now curdist is distance between init_plp.coords and polyline_coords[i-1]
			if(curdist >= distance)
				break; // polyline_coords[i-1] is far away enough, while polyline_coords[i] was still too close to init_coords
		}
		if(i == 0) {
			reached_polyline_extreme = true;
			return PolyLineGraph2D::polyline::pl_point(0,polyline_coords[0]);
		} else {
			ratio = (distance - prevdist) / (curdist - prevdist);
			return PolyLineGraph2D::polyline::pl_point(i-1,first_plus_ratio_of_segment(polyline_coords[i],polyline_coords[i-1], ratio));
		}
	} else if(direction == end) {
		if(init_plp.segment_index >= (polyline_coords.size() - 1)) {
			// return end
			reached_polyline_extreme = true;
			return PolyLineGraph2D::polyline::pl_point(polyline_coords.size() - 2,polyline_coords[polyline_coords.size() - 1]);
		}


		curdist = compute_2d_distance(polyline_coords[init_plp.segment_index+1],init_plp.coords);
		if(curdist >= distance) {
			ratio = distance / curdist;
			return PolyLineGraph2D::polyline::pl_point(init_plp.segment_index,first_plus_ratio_of_segment(init_plp.coords,polyline_coords[init_plp.segment_index+1], ratio));
		}
		for(i=init_plp.segment_index+1; i < polyline_coords.size() - 1; i++) {
			// segment polyline_coords[i+1] <--> polyline_coords[i]
			prevdist = curdist; // Now prevdist is distance between init_plp.coords and polyline_coords[i]
			curdist = compute_2d_distance(polyline_coords[i+1],init_plp.coords); // Now curdist is distance between init_plp.coords and polyline_coords[i+1]
			if(curdist >= distance)
				break; // polyline_coords[i+1] is far away enough, while polyline_coords[i] was still too close to init_coords
		}
		if(i == polyline_coords.size() - 1) {
			reached_polyline_extreme = true;
			return PolyLineGraph2D::polyline::pl_point(polyline_coords.size() - 2,polyline_coords[polyline_coords.size() - 1]);
		} else {
			ratio = (distance - prevdist) / (curdist - prevdist);
			return PolyLineGraph2D::polyline::pl_point(i,first_plus_ratio_of_segment(polyline_coords[i],polyline_coords[i+1], ratio));
		}
	} else
		std::invalid_argument("No start or end given as direction");
}

PolyLineGraph2D::polyline::pl_point PolyLineGraph2D::polyline::next_pl_point_by_length(const PolyLineGraph2D::polyline::pl_point init_plp, const ulong direction, const float length, bool &reached_polyline_extreme) {
	float prevlen, curlen;
	reached_polyline_extreme = false;
	float ratio;
	ulong i;

	if(direction == start) {
		curlen = compute_2d_distance(polyline_coords[init_plp.segment_index],init_plp.coords);
		if(curlen >= length) {
			ratio = length / curlen;
			return PolyLineGraph2D::polyline::pl_point(init_plp.segment_index,first_plus_ratio_of_segment(init_plp.coords,polyline_coords[init_plp.segment_index], ratio));
		}
		for(i=init_plp.segment_index; i > 0; i--) {
			// segment polyline_coords[i-1] <--> polyline_coords[i]
			prevlen = curlen; // Now prevdist is distance between init_plp.coords and polyline_coords[i]
			curlen += compute_2d_distance(polyline_coords[i-1],init_plp.coords); // Now curdist is distance between init_plp.coords and polyline_coords[i-1]
			if(curlen >= length)
				break; // polyline_coords[i-1] is far away enough, while polyline_coords[i] was still too close to init_coords
		}
		if(i == 0) {
			reached_polyline_extreme = true;
			return PolyLineGraph2D::polyline::pl_point(0,polyline_coords[0]);
		} else {
			ratio = (length - prevlen) / (curlen - prevlen);
			return PolyLineGraph2D::polyline::pl_point(i-1,first_plus_ratio_of_segment(polyline_coords[i],polyline_coords[i-1], ratio));
		}
	} else if(direction == end) {
		curlen = compute_2d_distance(polyline_coords[init_plp.segment_index+1],init_plp.coords);
		if(curlen >= length) {
			ratio = length / curlen;
			return PolyLineGraph2D::polyline::pl_point(init_plp.segment_index,first_plus_ratio_of_segment(init_plp.coords,polyline_coords[init_plp.segment_index+1], ratio));
		}
		for(i=init_plp.segment_index+1; i < polyline_coords.size() - 1; i++) {
			// segment polyline_coords[i+1] <--> polyline_coords[i]
			prevlen = curlen; // Now prevdist is distance between init_plp.coords and polyline_coords[i]
			curlen += compute_2d_distance(polyline_coords[i+1],init_plp.coords); // Now curdist is distance between init_plp.coords and polyline_coords[i+1]
			if(curlen >= length)
				break; // polyline_coords[i+1] is far away enough, while polyline_coords[i] was still too close to init_coords
		}
		if(i == polyline_coords.size() - 1) {
			reached_polyline_extreme = true;
			return PolyLineGraph2D::polyline::pl_point(polyline_coords.size() - 2,polyline_coords[polyline_coords.size() - 1]);
		} else {
			ratio = (length - prevlen) / (curlen - prevlen);
			return PolyLineGraph2D::polyline::pl_point(i,first_plus_ratio_of_segment(polyline_coords[i],polyline_coords[i+1], ratio));
		}
	} else
		std::invalid_argument("No start or end given as direction");
}
PolyLineGraph2D::polyline::pl_point PolyLineGraph2D::polyline::next_pl_point_by_distance(const ulong starting_extreme, const glm::vec2 coords, const float distance, bool &reached_polyline_extreme) {
	float prevdist, curdist;
	reached_polyline_extreme = false;
	float ratio;

	ulong i;

	ulong direction = get_other_end(starting_extreme);


	if(direction == start) {
		i = polyline_coords.size()-1;
		curdist = compute_2d_distance(polyline_coords[i],coords);

		if(curdist >= distance)
			return PolyLineGraph2D::polyline::pl_point(i-1,polyline_coords[i]);

		for(; i > 0; i--) {
			// segment polyline_coords[i-1] <--> polyline_coords[i]
			prevdist = curdist; // Now prevdist is distance between coords and polyline_coords[i]
			curdist = compute_2d_distance(polyline_coords[i-1],coords); // Now curdist is distance between coords and polyline_coords[i-1]
			if(curdist >= distance)
				break; // polyline_coords[i-1] is far away enough, while polyline_coords[i] was still too close to init_coords
		}

		if(i == 0) {
			reached_polyline_extreme = true;
			return PolyLineGraph2D::polyline::pl_point(0,polyline_coords[0]);
		} else {
			ratio = (distance - prevdist) / (curdist - prevdist);
			return PolyLineGraph2D::polyline::pl_point(i-1,first_plus_ratio_of_segment(polyline_coords[i],polyline_coords[i-1], ratio));
		}
	} else if(direction == end) {
		i = 0;
		curdist = compute_2d_distance(polyline_coords[i],coords);

		if(curdist >= distance)
			return PolyLineGraph2D::polyline::pl_point(i,polyline_coords[i]);

		for(; i < polyline_coords.size() - 1; i++) {
			// segment polyline_coords[i+1] <--> polyline_coords[i]
			prevdist = curdist; // Now prevdist is distance between init_plp.coords and polyline_coords[i]
			curdist = compute_2d_distance(polyline_coords[i+1],coords); // Now curdist is distance between init_plp.coords and polyline_coords[i+1]
			if(curdist >= distance)
				break; // polyline_coords[i+1] is far away enough, while polyline_coords[i] was still too close to init_coords
		}
		if(i == polyline_coords.size() - 1) {
			reached_polyline_extreme = true;
			return PolyLineGraph2D::polyline::pl_point(polyline_coords.size() - 2,polyline_coords[polyline_coords.size() - 1]);
		} else {
			ratio = (distance - prevdist) / (curdist - prevdist);
			return PolyLineGraph2D::polyline::pl_point(i,first_plus_ratio_of_segment(polyline_coords[i],polyline_coords[i+1], ratio));
		}
	} else
		std::invalid_argument("No start or end given as direction");
}

vector<PolyLineGraph2D::polyline::pl_point> PolyLineGraph2D::polyline::next_pl_points_by_distance(const PolyLineGraph2D::polyline::pl_point init_plp, const ulong direction, const float distance) {
	vector<PolyLineGraph2D::polyline::pl_point> res;
	pl_point next = init_plp;
	bool reached_polyline_extreme = false;

	while(!reached_polyline_extreme) {
		next = next_pl_point_by_distance(next,direction,distance,reached_polyline_extreme);
		res.push_back(pl_point(next));
	}

	return res;
}

vector<PolyLineGraph2D::polyline::pl_point> PolyLineGraph2D::polyline::split_equal_size_intervals(const ulong starting_extreme, const float distance) {
	ulong direction = get_other_end(starting_extreme);
	pl_point starting_plp = get_extreme_plp(starting_extreme);
	vector<PolyLineGraph2D::polyline::pl_point> t_res = next_pl_points_by_distance(starting_plp, direction, distance);
	vector<PolyLineGraph2D::polyline::pl_point> res;
	res.push_back(starting_plp);
	for(auto &plp:t_res)
		res.push_back(plp);
	return res;
}

void PolyLineGraph2D::polyline::next_pl_point_by_line_intersection(const pl_point init_plp, const ulong direction, const glm::vec3 &line, PolyLineGraph2D::polyline::pl_point &next, bool &found_quasiparallel_segment, PolyLineGraph2D::polyline::pl_point &next_before_quasiparallel_segment, bool &reached_polyline_extreme, bool &found) {
	ulong i;
	glm::vec4 cur_segm;
	glm::vec2 intersection;
	bool parallel,overlapped,intersection_found,quasiparallel_within_distance,valid;
	float distance;

	found_quasiparallel_segment = false;
	reached_polyline_extreme = false;
	found = false;

	if(direction == start) {
		cur_segm = glm::vec4(init_plp.coords,polyline_coords[init_plp.segment_index]);
		intersect_segment_line_no_quasiparallel(cur_segm, line, POLYLINE_NEXT_BY_LINE_INTERSECTION_MAX_QUASIPARALLEL_ANGLE_COS, POLYLINE_NEXT_BY_LINE_INTERSECTION_MAX_QUASIPARALLEL_DIST, parallel, overlapped, intersection_found, quasiparallel_within_distance, valid, distance, intersection);
		if(quasiparallel_within_distance) {
			next_before_quasiparallel_segment = init_plp;
			found_quasiparallel_segment = true;
			found = false;
			return;
		} else if (intersection_found) {
			next = PolyLineGraph2D::polyline::pl_point(init_plp.segment_index,intersection);
			found = true;
			return;
		}

		for(i=init_plp.segment_index; i > 0; i--) {
			// segment polyline_coords[i-1] <--> polyline_coords[i]
			cur_segm = glm::vec4(polyline_coords[i],polyline_coords[i-1]);
			intersect_segment_line_no_quasiparallel(cur_segm, line, POLYLINE_NEXT_BY_LINE_INTERSECTION_MAX_QUASIPARALLEL_ANGLE_COS, POLYLINE_NEXT_BY_LINE_INTERSECTION_MAX_QUASIPARALLEL_DIST, parallel, overlapped, intersection_found, quasiparallel_within_distance, valid, distance, intersection);

			if(quasiparallel_within_distance) {
				next_before_quasiparallel_segment = PolyLineGraph2D::polyline::pl_point(i-1, polyline_coords[i]);
				found_quasiparallel_segment = true;
				found = false;
				return;
			} else if (intersection_found) {
				next = PolyLineGraph2D::polyline::pl_point(i-1,intersection);
				found = true;
				return;
			}
		}

		// i = 0, reached start
		reached_polyline_extreme = true;
		found = false;
		return;

	} else if(direction == end) {

		cur_segm = glm::vec4(init_plp.coords,polyline_coords[init_plp.segment_index+1]);
		intersect_segment_line_no_quasiparallel(cur_segm, line, POLYLINE_NEXT_BY_LINE_INTERSECTION_MAX_QUASIPARALLEL_ANGLE_COS, POLYLINE_NEXT_BY_LINE_INTERSECTION_MAX_QUASIPARALLEL_DIST, parallel, overlapped, intersection_found, quasiparallel_within_distance, valid, distance, intersection);
		if(quasiparallel_within_distance) {
			next_before_quasiparallel_segment = init_plp;
			found_quasiparallel_segment = true;
			found = false;
			return;
		} else if (intersection_found) {
			next = PolyLineGraph2D::polyline::pl_point(init_plp.segment_index,intersection);
			found = true;
			return;
		}

		for(i=init_plp.segment_index+1; i < polyline_coords.size() - 1; i++) {
			// segment polyline_coords[i+1] <--> polyline_coords[i]
			cur_segm = glm::vec4(polyline_coords[i],polyline_coords[i+1]);
			intersect_segment_line_no_quasiparallel(cur_segm, line, POLYLINE_NEXT_BY_LINE_INTERSECTION_MAX_QUASIPARALLEL_ANGLE_COS, POLYLINE_NEXT_BY_LINE_INTERSECTION_MAX_QUASIPARALLEL_DIST, parallel, overlapped, intersection_found, quasiparallel_within_distance, valid, distance, intersection);

			if(quasiparallel_within_distance) {
				next_before_quasiparallel_segment = PolyLineGraph2D::polyline::pl_point(i, polyline_coords[i]);
				found_quasiparallel_segment = true;
				found = false;
				return;
			} else if (intersection_found) {
				next = PolyLineGraph2D::polyline::pl_point(i,intersection);
				found = true;
				return;
			}
		}
		// i = polyline_coords.size() - 1, reached end
		reached_polyline_extreme = true;
		found = false;
		return;

	} else
		std::invalid_argument("No start or end given as direction");
}

void PolyLineGraph2D::polyline::next_pl_point_by_line_intersection_bounded_distance(const pl_point init_plp, const ulong direction, const glm::vec3 &line, const float min_dist, const float max_dist, PolyLineGraph2D::polyline::pl_point &next, bool &found_quasiparallel_segment, PolyLineGraph2D::polyline::pl_point &next_before_quasiparallel_segment, bool &reached_polyline_extreme, bool &bounded_distance_violated, bool &found) {
	ulong i;
	glm::vec4 cur_segm;
	glm::vec2 intersection;
	bool parallel,overlapped,intersection_found,quasiparallel_within_distance,valid;
	float distance;

	bounded_distance_violated = false;
	found_quasiparallel_segment = false;
	reached_polyline_extreme = false;
	found = false;

	if(direction == start) {
		cur_segm = glm::vec4(init_plp.coords,polyline_coords[init_plp.segment_index]);
		intersect_segment_line_no_quasiparallel(cur_segm, line, POLYLINE_NEXT_BY_LINE_INTERSECTION_MAX_QUASIPARALLEL_ANGLE_COS, POLYLINE_NEXT_BY_LINE_INTERSECTION_MAX_QUASIPARALLEL_DIST, parallel, overlapped, intersection_found, quasiparallel_within_distance, valid, distance, intersection);
		if(quasiparallel_within_distance) {
			next_before_quasiparallel_segment = init_plp;
			found_quasiparallel_segment = true;
			found = false;
			return;
		} else if (intersection_found) {
			next = PolyLineGraph2D::polyline::pl_point(init_plp.segment_index,intersection);
			found = true;

			float dsq = squared_2d_distance(intersection,init_plp.coords);
			if(dsq < (min_dist*min_dist) || dsq > (max_dist*max_dist)) {
				bounded_distance_violated = true;
				found = false;
			}

			return;
		}

		for(i=init_plp.segment_index; i > 0; i--) {
			// segment polyline_coords[i-1] <--> polyline_coords[i]
			cur_segm = glm::vec4(polyline_coords[i],polyline_coords[i-1]);
			intersect_segment_line_no_quasiparallel(cur_segm, line, POLYLINE_NEXT_BY_LINE_INTERSECTION_MAX_QUASIPARALLEL_ANGLE_COS, POLYLINE_NEXT_BY_LINE_INTERSECTION_MAX_QUASIPARALLEL_DIST, parallel, overlapped, intersection_found, quasiparallel_within_distance, valid, distance, intersection);

			if(quasiparallel_within_distance) {
				next_before_quasiparallel_segment = PolyLineGraph2D::polyline::pl_point(i-1, polyline_coords[i]);
				found_quasiparallel_segment = true;
				found = false;
				return;
			} else if (intersection_found) {
				next = PolyLineGraph2D::polyline::pl_point(i-1,intersection);
				found = true;

				float dsq = squared_2d_distance(intersection,init_plp.coords);
				if(dsq < (min_dist*min_dist) || dsq > (max_dist*max_dist)) {
					bounded_distance_violated = true;
					found = false;
				}

				return;
			}
		}

		// i = 0, reached start
		reached_polyline_extreme = true;
		found = false;
		return;

	} else if(direction == end) {

		cur_segm = glm::vec4(init_plp.coords,polyline_coords[init_plp.segment_index+1]);
		intersect_segment_line_no_quasiparallel(cur_segm, line, POLYLINE_NEXT_BY_LINE_INTERSECTION_MAX_QUASIPARALLEL_ANGLE_COS, POLYLINE_NEXT_BY_LINE_INTERSECTION_MAX_QUASIPARALLEL_DIST, parallel, overlapped, intersection_found, quasiparallel_within_distance, valid, distance, intersection);
		if(quasiparallel_within_distance) {
			next_before_quasiparallel_segment = init_plp;
			found_quasiparallel_segment = true;
			found = false;
			return;
		} else if (intersection_found) {
			next = PolyLineGraph2D::polyline::pl_point(init_plp.segment_index,intersection);
			found = true;

			float dsq = squared_2d_distance(intersection,init_plp.coords);
			if(dsq < (min_dist*min_dist) || dsq > (max_dist*max_dist)) {
				bounded_distance_violated = true;
				found = false;
			}

			return;
		}

		for(i=init_plp.segment_index+1; i < polyline_coords.size() - 1; i++) {
			// segment polyline_coords[i+1] <--> polyline_coords[i]
			cur_segm = glm::vec4(polyline_coords[i],polyline_coords[i+1]);
			intersect_segment_line_no_quasiparallel(cur_segm, line, POLYLINE_NEXT_BY_LINE_INTERSECTION_MAX_QUASIPARALLEL_ANGLE_COS, POLYLINE_NEXT_BY_LINE_INTERSECTION_MAX_QUASIPARALLEL_DIST, parallel, overlapped, intersection_found, quasiparallel_within_distance, valid, distance, intersection);

			if(quasiparallel_within_distance) {
				next_before_quasiparallel_segment = PolyLineGraph2D::polyline::pl_point(i, polyline_coords[i]);
				found_quasiparallel_segment = true;
				found = false;
				return;
			} else if (intersection_found) {
				next = PolyLineGraph2D::polyline::pl_point(i,intersection);
				found = true;

				float dsq = squared_2d_distance(intersection,init_plp.coords);
				if(dsq < (min_dist*min_dist) || dsq > (max_dist*max_dist)) {
					bounded_distance_violated = true;
					found = false;
				}

				return;
			}
		}
		// i = polyline_coords.size() - 1, reached end
		reached_polyline_extreme = true;
		found = false;
		return;

	} else
		std::invalid_argument("No start or end given as direction");
}

ulong PolyLineGraph2D::polyline::get_amount_of_segments() {
	return polyline_coords.size()-1;
}

PolyLineGraph2D::polyline::polyline() {}

PolyLineGraph2D::polyline::polyline(const ulong s, const ulong e, const vector<glm::vec2> &pcs) : start(s), end(e), polyline_coords(pcs) {
	update_length();
}

PolyLineGraph2D::polyline::polyline(const ulong &s, const ulong &e, const vector<glm::vec2> &pcs, float length) : start(s), end(e), polyline_coords(pcs), length(length) { }

PolyLineGraph2D::polyline::~polyline() {
	//delete polyline_coords;
}

#define SQROOT_TWO 1.414
#define PL_CELL_SPLIT_RATIO (SQROOT_TWO+0.1)

vector<pair<ulong,ulong>> PolyLineGraph2D::polyline::get_intersectedcells_2dmap_vec(const float cell_dim) {
	vector<pair<ulong,ulong>> intersectedcells;
	vector<PolyLineGraph2D::polyline::pl_point> plps = split_equal_size_intervals(start,cell_dim / PL_CELL_SPLIT_RATIO);
	pair<ulong,ulong> prev_cell;
	pair<ulong,ulong> cur_cell;
	bool on_boundary;

	for(const auto &plp: plps) {
		cur_cell = get_2dmap_cell_from_coords(cell_dim, plp.coords, on_boundary);
		if(!on_boundary && (intersectedcells.size() == 0 || prev_cell != cur_cell)) {
				intersectedcells.push_back(cur_cell);
				prev_cell = cur_cell;
			}
	}

	return intersectedcells;
}

set<pair<ulong,ulong>> PolyLineGraph2D::polyline::get_intersectedcells_2dmap_set(const float cell_dim) {
	set<pair<ulong,ulong>> intersectedcells;
	vector<PolyLineGraph2D::polyline::pl_point> plps = split_equal_size_intervals(start,cell_dim / PL_CELL_SPLIT_RATIO);
	pair<ulong,ulong> prev_cell;
	pair<ulong,ulong> cur_cell;
	bool on_boundary;

	for(const auto &plp: plps) {
		cur_cell = get_2dmap_cell_from_coords(cell_dim, plp.coords, on_boundary);
		if(!on_boundary && (intersectedcells.size() == 0 || prev_cell != cur_cell)) {
				intersectedcells.insert(cur_cell);
				prev_cell = cur_cell;
			}
	}

	return intersectedcells;
}

bool PolyLineGraph2D::polyline::connects(const ulong a) {
	return a == start || a == end;
}

bool PolyLineGraph2D::polyline::connects(const ulong a, const ulong b) {
	return (a == start && b == end) || (b == start && a == end);
}

float PolyLineGraph2D::polyline::compute_distancesq(const glm::vec2 &p_coords, ulong &closest_segm, glm::vec2 &projection) {
	float min_dist = minimum_distancesq(p_coords,polyline_coords[0],polyline_coords[1],projection);
	ulong min_segm = 0;
	closest_segm = 0;
	float cur_dist;
	glm::vec2 cur_projection;

	for(ulong i=2; i < polyline_coords.size(); i++) {
		cur_dist = minimum_distancesq(p_coords,polyline_coords[i-1],polyline_coords[i],cur_projection);
		if(cur_dist < min_dist) {
			min_dist = cur_dist;
			projection = cur_projection;
			closest_segm = i-1;
		}
	}

	return min_dist;
}

float PolyLineGraph2D::polyline::compute_distancesq(const glm::vec2 &p_coords, pl_point &plp) {
	return compute_distancesq(p_coords,plp.segment_index,plp.coords);
}

PolyLineGraph2D::plg_point::plg_point() {}

PolyLineGraph2D::plg_point::plg_point(const ulong input_polyline_id, const ulong input_segment_index, const glm::vec2 input_coords) :
		polyline_id(input_polyline_id), plp(polyline::pl_point(input_segment_index,input_coords)) {}

PolyLineGraph2D::plg_point::plg_point(const ulong input_polyline_id, const polyline::pl_point &plp) :
		polyline_id(input_polyline_id), plp(plp) {}
PolyLineGraph2D::plg_point::plg_point(const pair<ulong,ulong> &input_polyline_segment_index, const glm::vec2 input_coords) :
				polyline_id(input_polyline_segment_index.first), plp(polyline::pl_point(input_polyline_segment_index.second,input_coords)) {}

bool PolyLineGraph2D::plg_point::operator==(const plg_point& a) const {
	return polyline_id == a.polyline_id && plp == a.plp;
}

bool PolyLineGraph2D::plgs_point::operator==(const plgs_point& a) const {
	return plg_id == a.plg_id && plgp == a.plgp;
}

PolyLineGraph2D::plgs_point::plgs_point() {}

PolyLineGraph2D::plgs_point::plgs_point(const ulong plg_id, const plg_point &plgp) :
		plg_id(plg_id), plgp(plgp) {}

vector<glm::vec4> PolyLineGraph2D::polyline::get_segments_list() {
	vector<glm::vec4> res;

	for(ulong i = 1; i < polyline_coords.size(); i++) {
		res.push_back(glm::vec4(polyline_coords[i-1],polyline_coords[i]));
	}

	return res;
}

ulong PolyLineGraph2D::polyline::get_other_end(ulong extreme) {
	if(extreme == start)
		return end;
	else if(extreme == end)
		return start;
	else
		std::invalid_argument("No start or end given as direction");
}

bool linearizable_polyline(const vector<glm::vec2> &polyline_coords, const ulong start, const ulong end, const float max_linearizability_distsq) {
	glm::vec3 line = compute_2dline(polyline_coords[start],polyline_coords[end]);

	for(ulong i = start+1; i < end; i++)
		if(distance_point_line_sq(polyline_coords[i], line) > max_linearizability_distsq)
			// current ppline is not compatible with computed 2D line
			return false;

	return true;
}

ulong find_max_se(const vector<glm::vec2> &polyline_coords, const ulong start, const ulong max_se, const float max_linearizability_distsq) {
	if(max_se <= start)
		return start;
	for(ulong cur_se = max_se; cur_se > start+1; cur_se--)
		if(linearizable_polyline(polyline_coords,start,cur_se,max_linearizability_distsq))
			return cur_se;
	return start+1;
}

ulong find_min_eb(const vector<glm::vec2> &polyline_coords, const ulong end, const ulong min_eb, const float max_linearizability_distsq) {
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
pair<ulong,ulong> find_compatible_se_eb(const vector<glm::vec2> &polyline_coords, const ulong start, const ulong end, const float max_linearizability_distsq) {
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

vector<glm::vec2> simplify_polyline(const vector<glm::vec2> &polyline_coords, const float max_linearizability_dist) {
	ulong start,end;
	pair<ulong,ulong> se_eb;
	ulong se, eb;
	ulong max_se, min_eb;
	vector<glm::vec2> simplified_polyline_coords, simplified_polyline_coords_end;

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

	for(vector<glm::vec2>::reverse_iterator it = simplified_polyline_coords_end.rbegin(); it != simplified_polyline_coords_end.rend(); it++)
		simplified_polyline_coords.push_back(*it);

	return simplified_polyline_coords;
}

void PolyLineGraph2D::polyline::simplify() {
	simplify(MAXIMUM_LINEARIZABILITY_DISTANCE);
}

void PolyLineGraph2D::polyline::simplify(const float max_linearizability_dist) {
	vector<glm::vec2> spl = simplify_polyline(polyline_coords,max_linearizability_dist);
	polyline_coords.clear();
	polyline_coords = spl;
	update_length();
}

bool PolyLineGraph2D::polyline::operator==(const polyline& p) const
{
    return (start == p.start && end == p.end && vec_glm_vec2_equal(polyline_coords,p.polyline_coords)) || (start == p.end && end == p.start && vec_glm_vec2_equal_inv(polyline_coords,p.polyline_coords));
}

bool PolyLineGraph2D::polyline::is_loop() {
	return start == end;
}

void PolyLineGraph2D::polyline::invalidate() {
	clear_coords();
	length = INVALID_POLYLINE_LENGTH;
}

void PolyLineGraph2D::remove_connection(const ulong node_id, const ulong polyline_id) {
	connections[node_id].erase( std::remove( connections[node_id].begin(), connections[node_id].end(), polyline_id ), connections[node_id].end() );

	if(connections[node_id].size() == 0)
		invalidate_node(node_id);
}

void PolyLineGraph2D::remove_polyline(ulong polyline_id) {
	PolyLineGraph2D::polyline &p = polylines[polyline_id];
	remove_connection(p.start,polyline_id);
	remove_connection(p.end,polyline_id);
	p.invalidate();
}

void PolyLineGraph2D::polyline::clear_coords() {
	polyline_coords.clear();
	polyline_coords.resize(0);
	update_length();
}

PolyLineGraph2D::polyline PolyLineGraph2D::polyline::merge_polylines(const polyline& p1,const polyline& p2) {
	ulong common_extreme;
	ulong p3s,p3e;
	vector<glm::vec2> p3coords;

	if(p1.start == p2.start) {
		//  <---- P1 P2 ---->
		p3s = p1.end;
		p3e = p2.end;
		for(std::vector<glm::vec2>::const_reverse_iterator rit = p1.polyline_coords.rbegin(); rit!= p1.polyline_coords.rend(); ++rit)
			p3coords.push_back(*rit);
		std::vector<glm::vec2>::const_iterator it = p2.polyline_coords.begin();
		it++;
		for (; it != p2.polyline_coords.end(); ++it)
			p3coords.push_back(*it);
		return polyline(p3s,p3e,p3coords);
	} else if(p1.start == p2.end) {
		//  P2 ----> P1 ---->

		p3s = p2.start;
		p3e = p1.end;
		for(std::vector<glm::vec2>::const_iterator it = p2.polyline_coords.begin(); it!= p2.polyline_coords.end(); ++it)
			p3coords.push_back(*it);
		std::vector<glm::vec2>::const_iterator it = p1.polyline_coords.begin();
		it++;
		for (; it != p1.polyline_coords.end(); ++it)
			p3coords.push_back(*it);
		return polyline(p3s,p3e,p3coords);

	} else if(p1.end == p2.start) {
		//  P1 ----> P2 ---->

		p3s = p1.start;
		p3e = p2.end;
		for(std::vector<glm::vec2>::const_iterator it = p1.polyline_coords.begin(); it!= p1.polyline_coords.end(); ++it)
			p3coords.push_back(*it);
		std::vector<glm::vec2>::const_iterator it = p2.polyline_coords.begin();
		it++;
		for (; it != p2.polyline_coords.end(); ++it)
			p3coords.push_back(*it);
		return polyline(p3s,p3e,p3coords);

	} else if(p1.end == p2.end) {
		//  P1 ----> <---- P2
		p3s = p1.start;
		p3e = p2.start;
		for(std::vector<glm::vec2>::const_iterator it = p1.polyline_coords.begin(); it!= p1.polyline_coords.end(); ++it)
			p3coords.push_back(*it);
		std::vector<glm::vec2>::const_reverse_iterator rit = p2.polyline_coords.rbegin();
		rit++;
		for (; rit != p2.polyline_coords.rend(); ++rit)
			p3coords.push_back(*rit);
		return polyline(p3s,p3e,p3coords);
	} else
		throw std::invalid_argument( "cannot merge disconnected polylines" );
}

void PolyLineGraph2D::simplify()
{
	simplify(MAXIMUM_LINEARIZABILITY_DISTANCE);
}

void PolyLineGraph2D::simplify(const float max_linearizability_dist)
{
	for(ulong pl_id=0; pl_id < polylines.size(); pl_id++)
		if(is_valid_polyline(pl_id))
			polylines[pl_id].simplify(max_linearizability_dist);
}

void PolyLineGraph2D::invalidate_node(ulong node_id) {
	nodes_coords[node_id] = glm::vec2(INVALID_POINT_COORDS,INVALID_POINT_COORDS);
	for(const auto pid : connections[node_id])
		remove_polyline(pid);
	connections[node_id].clear();
	connections[node_id].resize(0);
}

bool PolyLineGraph2D::is_valid_node(const ulong node_id) {
	return nodes_coords[node_id][0] != INVALID_POINT_COORDS && nodes_coords[node_id][1] != INVALID_POINT_COORDS;
}

bool PolyLineGraph2D::is_valid_polyline(const ulong polyline_id) {
	return is_valid_node(polylines[polyline_id].start)
			&& is_valid_node(polylines[polyline_id].end)
			&& polylines[polyline_id].polyline_coords.size() > 1
			&& get_node_coords(polylines[polyline_id].start) == polylines[polyline_id].polyline_coords[0]
			&& get_node_coords(polylines[polyline_id].end) == polylines[polyline_id].polyline_coords[polylines[polyline_id].polyline_coords.size()-1];
}

void PolyLineGraph2D::optimize()
{
	const ulong polylines_amount = polylines.size();
	for(ulong i=0; i < polylines_amount; i++)
		if(is_valid_polyline(i))
			polylines[i].simplify();
}

vector<glm::vec4> PolyLineGraph2D::get_segments_list() {
	vector<glm::vec4> res;

	for(ulong i=0; i < polylines.size(); i++)
		if(is_valid_polyline(i)) {
			const polyline &p = polylines[i];
			const vector<glm::vec4> cur_res = p.get_segments_list();
			for(const auto &s : cur_res)
				res.push_back(s);
		}


	return res;
}

vector<vector<glm::vec4>> PolyLineGraph2D::get_segments_grouped_by_polyline() {
	vector<vector<glm::vec4>> res;

	for(ulong i=0; i < polylines.size(); i++)
		if(is_valid_polyline(i)) {
			const polyline &p = polylines[i];
			res.push_back(p.get_segments_list());
		}


	return res;
}

vector<pair<ulong,vector<glm::vec4>>> PolyLineGraph2D::get_segments_grouped_by_polyline_with_polyline_ids() {
	vector<pair<ulong,vector<glm::vec4>>> res;

	for(ulong i=0; i < polylines.size(); i++)
		if(is_valid_polyline(i)) {
			const polyline &p = polylines[i];
			res.push_back(make_pair(i,p.get_segments_list()));
		}

	return res;
}

vector<vector<glm::vec4>> PolyLineGraph2D::get_segments_grouped_by_component() {
	vector<vector<glm::vec4>> res;
	ulong other_end;
	vector<glm::vec4> cur_pl_res;

	bool* explored = new bool[get_nodes_amount()];
	for(ulong i=0; i < get_nodes_amount(); i++)
		explored[i] = false;

	bool* in_to_explore = new bool[get_nodes_amount()];
	for(ulong i=0; i < get_nodes_amount(); i++)
		in_to_explore[i] = false;

	stack<ulong> to_explore;

	// Explore
	for(ulong start_node = 0; start_node < get_nodes_amount(); start_node++) {
		if(!explored[start_node] && is_valid_node(start_node)) {
			vector<glm::vec4> cur_res;
			explored[start_node] = true;

			for(const auto p_id : connections[start_node]) {
				if(is_valid_polyline(p_id))
				{
					other_end = polylines[p_id].get_other_end(start_node);
					cur_pl_res = polylines[p_id].get_segments_list();
					for(const auto &s : cur_pl_res)
						cur_res.push_back(s);
					in_to_explore[other_end] = true;
					to_explore.push(other_end);
				}
			}

			while(!to_explore.empty()) {
				const ulong cur_node = to_explore.top();
				to_explore.pop();
				in_to_explore[cur_node] = false;
				explored[cur_node] = true;

				for(const auto p_id : connections[cur_node]) {
					if(is_valid_polyline(p_id))
					{
						other_end = polylines[p_id].get_other_end(cur_node);
						// " || other_end == cur_node " deals with 1-loops
						if(!explored[other_end] || other_end == cur_node) {
							cur_pl_res = polylines[p_id].get_segments_list();
							for(const auto &s : cur_pl_res)
								cur_res.push_back(s);
							if(!in_to_explore[other_end] && other_end != cur_node)
								to_explore.push(other_end);
					}
					}
				}
			}

			res.push_back(cur_res);
		}
	}

	return res;
}

glm::vec2 PolyLineGraph2D::get_node_coords(ulong node_id) {
	return nodes_coords[node_id];
}

vector<glm::vec2> PolyLineGraph2D::get_nodes_coords() {
	return nodes_coords;
}

ulong PolyLineGraph2D::get_nodes_amount() {
	return nodes_amount;
}

ulong PolyLineGraph2D::get_polylines_amount() {
	return polylines.size();
}

void PolyLineGraph2D::add_node_coords(const glm::vec2 &p_coords)
{
	nodes_coords.push_back(p_coords);
}

vector<PolyLineGraph2D::plg_point> convert_vec_pl_point_to_plg_point(const vector<PolyLineGraph2D::polyline::pl_point> &vpl, const ulong polyline_id) {
	vector<PolyLineGraph2D::plg_point> res;

	for(auto &pl: vpl)
		res.push_back(PolyLineGraph2D::plg_point(polyline_id,pl));

	return res;
}

std::ostream &operator<< (std::ostream &out, const PolyLineGraph2D::plg_point &plgp) {
	out << "PLGP(" << plgp.polyline_id << "," << plgp.plp.segment_index << "," << plgp.plp.coords << ")";
	return out;
}

new_3dpoint_and_sides_plgp_matches make_new_3dpoint_and_sides_plgp_matches(vector<new_3dpoint_plgp_matches> valid_points_direction1,vector<ulong> directions1,new_3dpoint_plgp_matches central_point,vector<new_3dpoint_plgp_matches> valid_points_direction2,vector<ulong> directions2) {
	return make_tuple(make_pair(valid_points_direction1,directions1),central_point,make_pair(valid_points_direction2,directions2));
}

vector<new_3dpoint_plgp_matches> new_3dpoint_and_sides_plgp_matches_to_vector(const new_3dpoint_and_sides_plgp_matches &p3d_with_sides) {
	vector<new_3dpoint_plgp_matches> res;
	for(int i=get<0>(p3d_with_sides).first.size()-1; i >= 0; i--)
		res.push_back((get<0>(p3d_with_sides).first)[i]);
	res.push_back(get<1>(p3d_with_sides));
	for(int i=0; i < get<2>(p3d_with_sides).first.size(); i++)
		res.push_back((get<2>(p3d_with_sides).first)[i]);
	return res;
}

vector<glm::vec2> convert_vec_pl_point_to_vec2(const vector<PolyLineGraph2D::polyline::pl_point> &vpl) {
	vector<glm::vec2> res;
	for(const auto &pl : vpl)
		res.push_back(pl.coords);
	return res;
}

vector<pair<ulong,ulong>> find_closest_pairs_with_max_dist(pair<vector<ulong>,vector<glm::vec2>> nodes_data, const float max_dist) {
	const vector<ulong> &ids = nodes_data.first;
	const vector<glm::vec2> &p2ds = nodes_data.second;

	const float max_dist_sq = max_dist * max_dist;
	vector<ulong> closest_pt(p2ds.size());
	vector<pair<ulong,ulong>> res;
	float mindistsq,curdistsq;
	ulong min_id;

	for(ulong i=0;i<p2ds.size();i++) {
		mindistsq =  std::numeric_limits<float>::max();
		min_id=-1;

		for(ulong j=0;j<p2ds.size();j++)
			if (j!=i) {
				curdistsq = squared_2d_distance(p2ds[i],p2ds[j]);
				if(curdistsq < mindistsq) {
					mindistsq = curdistsq;
					min_id = j;
				}
			}

		closest_pt[i] = min_id;

		if(min_id < i) {
			// Than the closest point has been already computed for min_id

			if(i == closest_pt[min_id] && squared_2d_distance(p2ds[i],p2ds[min_id]) <= max_dist_sq)
				// Closest (reciprocal) pair
				res.push_back(make_pair(ids[i],ids[min_id]));
		}
	}

	closest_pt.clear();

	return res;
}

/**
 * Get close pair of extremes linked by a segment that follows the direction from both extremes
 */
vector<pair<ulong,ulong>> find_closest_pairs_with_max_dist_following_direction(tuple<vector<ulong>,vector<glm::vec2>,vector<glm::vec2>> nodes_data, const float max_dist, const float min_cos) {
	const vector<ulong> &ids = get<0>(nodes_data);
	const vector<glm::vec2> &p2ds = get<1>(nodes_data);
	const vector<glm::vec2> &dirs = get<2>(nodes_data);

	const float max_dist_sq = max_dist * max_dist;
	vector<ulong> closest_pt(p2ds.size());
	vector<pair<ulong,ulong>> res;
	float mindistsq,curdistsq;
	ulong min_id;

	for(ulong i=0;i<p2ds.size();i++) {
		mindistsq =  std::numeric_limits<float>::max();
		min_id=i;

		for(ulong j=0;j<p2ds.size();j++)
			if (j!=i) {
				curdistsq = squared_2d_distance(p2ds[i],p2ds[j]);
				if(curdistsq < mindistsq && curdistsq < max_dist_sq) {
					// compute_anglecos_vec2_vec2(polyline_coords[i-1],polyline_coords[i],polyline_coords[i-2],polyline_coords[i-1])
					//cout << "cos(" << p2ds[j]-p2ds[i] << "," << dirs[i] << ") = " << compute_anglecos_vec2_vec2(p2ds[j]-p2ds[i],dirs[i]) << "\n";

					if(abs(compute_anglecos_vec2_vec2(p2ds[j]-p2ds[i],dirs[i])) < min_cos)
						continue;
					if(abs(compute_anglecos_vec2_vec2(p2ds[j]-p2ds[i],dirs[j])) < min_cos)
						continue;

					mindistsq = curdistsq;
					min_id = j;
				}
			}

		closest_pt[i] = min_id;

		if(min_id < i) {
			// Than the closest point has been already computed for min_id

			if(i == closest_pt[min_id] && squared_2d_distance(p2ds[i],p2ds[min_id]) <= max_dist_sq)
				// Closest (reciprocal) pair
				res.push_back(make_pair(ids[i],ids[min_id]));
		}
	}

	closest_pt.clear();

	//cout << "Found pairs: " << res.size() << endl;

	return res;
}

vector<vector<ulong>> PolyLineGraph2D::get_connections() {
	return connections;
}

vector<PolyLineGraph2D::polyline> PolyLineGraph2D::get_polylines() {
	return polylines;
}

ulong PolyLineGraph2D::get_next_node_id() {
	ulong cur = next_node_id;
	next_node_id++;
	nodes_amount++;
	visited_nodes.push_back(false);
	return cur;
}

bool PolyLineGraph2D::has_loop(ulong node_id)
{
	for(const auto connection : connections[node_id])
		if(polylines[connection].is_loop())
			return true;
	return false;
}

/* Hub, either:
 * - >= 3 connected polylines
 * - 2 connected polylines, provided 1 polyline is a loop
 */
bool PolyLineGraph2D::is_hub(ulong node_id)
{
	const int connections_amount = connections[node_id].size();
	return (connections_amount > 2) || (connections_amount == 2 && has_loop(node_id));
}

/* Extreme, either:
 * - 1 connected polyline (not a loop)
 */
bool PolyLineGraph2D::is_extreme(ulong node_id)
{
	return connections[node_id].size() == 1 && !polylines[connections[node_id][0]].is_loop();
}

/* LoopNode, either:
 * - 1 connected polyline (a loop)
 */
bool PolyLineGraph2D::is_loopnode(ulong node_id)
{
	return connections[node_id].size() == 1 && polylines[connections[node_id][0]].is_loop();
}

vector<ulong> PolyLineGraph2D::get_hub_nodes()
{
	vector<ulong> res;
	for(ulong node_id = 0; node_id < get_nodes_amount(); node_id++)
		if(is_valid_node(node_id) && is_hub(node_id))
			res.push_back(node_id);
	return res;
}

vector<ulong> PolyLineGraph2D::get_extreme_nodes()
{
	vector<ulong> res;
	for(ulong node_id = 0; node_id < get_nodes_amount(); node_id++)
		if(is_valid_node(node_id) && is_extreme(node_id))
			res.push_back(node_id);
	return res;
}

pair<vector<ulong>,vector<glm::vec2>> PolyLineGraph2D::get_extreme_nodes_ids_and_coords()
{
	vector<ulong> res;
	vector<glm::vec2> res_coords;
	for(ulong node_id = 0; node_id < get_nodes_amount(); node_id++)
		if(is_valid_node(node_id) && is_extreme(node_id)) {
			res.push_back(node_id);
			res_coords.push_back(nodes_coords[node_id]);
		}

	return make_pair(res,res_coords);
}

tuple<vector<ulong>,vector<glm::vec2>,vector<glm::vec2>> PolyLineGraph2D::get_extreme_nodes_ids_and_coords_and_direction()
{
	vector<ulong> res;
	vector<glm::vec2> res_coords;
	vector<glm::vec2> res_dirs;
	for(ulong node_id = 0; node_id < get_nodes_amount(); node_id++)
		if(is_valid_node(node_id) && is_extreme(node_id)) {
			res.push_back(node_id);
			res_coords.push_back(nodes_coords[node_id]);
			res_dirs.push_back(polylines[connections[node_id][0]].get_extreme_direction_length_given_length(node_id,6).first);
		}

	return make_tuple(res,res_coords,res_dirs);
}

vector<ulong> PolyLineGraph2D::get_loopnodes()
{
	vector<ulong> res;
	for(ulong node_id = 0; node_id < get_nodes_amount(); node_id++)
		if(is_valid_node(node_id) && is_loopnode(node_id))
			res.push_back(node_id);
	return res;
}

vector<ulong> PolyLineGraph2D::get_nodes_with_loops()
{
	vector<ulong> res;
	for(ulong node_id = 0; node_id < get_nodes_amount(); node_id++)
		if(is_valid_node(node_id) && has_loop(node_id))
			res.push_back(node_id);
	return res;
}

vector<glm::vec2> PolyLineGraph2D::get_hub_nodes_coords()
{
	vector<glm::vec2> res;
	for(ulong node_id = 0; node_id < get_nodes_amount(); node_id++)
		if(is_valid_node(node_id) && is_hub(node_id))
			res.push_back(nodes_coords[node_id]);
	return res;
}

vector<glm::vec2> PolyLineGraph2D::get_extreme_nodes_coords()
{
	vector<glm::vec2> res;
	for(ulong node_id = 0; node_id < get_nodes_amount(); node_id++)
		if(is_valid_node(node_id) && is_extreme(node_id))
			res.push_back(nodes_coords[node_id]);
	return res;
}

vector<glm::vec2> PolyLineGraph2D::get_loopnodes_coords()
{
	vector<glm::vec2> res;
	for(ulong node_id = 0; node_id < get_nodes_amount(); node_id++)
		if(is_valid_node(node_id) && is_loopnode(node_id))
			res.push_back(nodes_coords[node_id]);
	return res;
}

vector<glm::vec2> PolyLineGraph2D::get_nodes_with_loops_coords()
{
	vector<glm::vec2> res;
	for(ulong node_id = 0; node_id < get_nodes_amount(); node_id++)
		if(is_valid_node(node_id) && has_loop(node_id))
			res.push_back(nodes_coords[node_id]);
	return res;
}

bool PolyLineGraph2D::is_connected_node(const ulong start_node,const ulong end_node) {
	vector<bool> visited_vec;
	visited_vec.push_back(start_node);
	visited_nodes[start_node]=true;

/*	set<ulong> visited;
	visited.insert(start_node);*/
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



PolyLineGraph2D::PolyLineGraph2D(const vector<polyline> &polylines, const vector<vector<ulong>> &connections, const vector<bool> &visited_nodes, const ulong &real_nodes_amount, const ulong &nodes_amount, const ulong &next_node_id, const vector<glm::vec2> &nodes_coords) :
		polylines(polylines), connections(connections), visited_nodes(visited_nodes), real_nodes_amount(real_nodes_amount), nodes_amount(nodes_amount), next_node_id(next_node_id), nodes_coords(nodes_coords) {}

PolyLineGraph2D::~PolyLineGraph2D() {}

void update_new_3dpoint_plgp_matches(new_3dpoint_plgp_matches &pt, const int new_view, const PolyLineGraph2D::plg_point &new_observation, const glm::vec3 &new_coords) {
	get<0>(pt) = new_coords;
	get<1>(pt).push_back(new_observation);
	get<2>(pt).push_back(new_view);
}

bool PolyLineGraph2D::is_connected_node(ulong start_node,ulong end_node, ulong max_jumps) {
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

bool can_follow_polyline_without_going_outside_radius(const vector<glm::vec2> &polyline_coords,const glm::vec2 &a,const float detection_radius_sq) {
	for(const auto &p : polyline_coords)
		if(squared_2d_distance(p,a) > detection_radius_sq)
			return false;

	return true;
}

bool PolyLineGraph2D::is_connected_node_inside_radius(const ulong start_node,const ulong end_node, const float detection_radius) {
	vector<bool> visited_vec;
	visited_vec.push_back(start_node);
	visited_nodes[start_node]=true;

	const glm::vec2 &start_coords = nodes_coords[start_node];

	const float detection_radius_sq = detection_radius*detection_radius;

	if(squared_2d_distance(nodes_coords[start_node],nodes_coords[end_node]) > detection_radius_sq)
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

bool PolyLineGraph2D::is_connected_polyline(ulong polyline_id_start,ulong polyline_id_end) {
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

bool PolyLineGraph2D::is_connected_polyline(ulong polyline_id_start,ulong polyline_id_end, ulong max_jumps) {
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
inline float min_dist_sq(const glm::vec2 &p,const glm::vec2 &x,const glm::vec2 &y) {
	return min(squared_2d_distance(p,x),squared_2d_distance(p,y));
}

bool can_follow_polyline_without_going_outside_radius(const vector<glm::vec2> &polyline_coords,const glm::vec2 &s_a,const glm::vec2 &s_b,const float detection_radius_sq) {
	for(const auto &p : polyline_coords)
		if(min_dist_sq(p,s_a,s_b) > detection_radius_sq)
			return false;

	return true;
}

bool PolyLineGraph2D::is_connected_polyline_inside_radius(ulong polyline_id_start, ulong polyline_id_end, const float detection_radius) {
	const polyline &p_start = polylines[polyline_id_start];
	const polyline &p_end = polylines[polyline_id_end];

	const float detection_radius_sq = detection_radius*detection_radius;

	const glm::vec2 &s_a = nodes_coords[p_start.start];
	const glm::vec2 &s_b = nodes_coords[p_start.end];
	const glm::vec2 &e_a = nodes_coords[p_end.start];
	const glm::vec2 &e_b = nodes_coords[p_end.end];

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

vector<ulong> PolyLineGraph2D::compute_components_node_ids() {
	vector<ulong> components_node_ids(get_nodes_amount());

	ulong other_end;

	bool* explored = new bool[get_nodes_amount()];
	for(ulong i=0; i < get_nodes_amount(); i++)
		explored[i] = false;

	bool* in_to_explore = new bool[get_nodes_amount()];
	for(ulong i=0; i < get_nodes_amount(); i++)
		in_to_explore[i] = false;

	stack<ulong> to_explore;

	ulong cur_component_id = 0;

	// Explore
	for(ulong start_node = 0; start_node < get_nodes_amount(); start_node++) {
		if(!explored[start_node]) {
			vector<glm::vec4> cur_res;
			explored[start_node] = true;

			components_node_ids[start_node] = cur_component_id;

			for(const auto p_id : connections[start_node]) {
				other_end = polylines[p_id].get_other_end(start_node);
				in_to_explore[other_end] = true;
				to_explore.push(other_end);
			}

			while(!to_explore.empty()) {
				const ulong cur_node = to_explore.top();
				to_explore.pop();

				components_node_ids[cur_node] = cur_component_id;

				in_to_explore[cur_node] = false;
				explored[cur_node] = true;

				for(const auto p_id : connections[cur_node]) {
					other_end = polylines[p_id].get_other_end(cur_node);
					// " || other_end == cur_node " deals with 1-loops
					if(!explored[other_end] || other_end == cur_node) {
						if(!in_to_explore[other_end] && other_end != cur_node)
							to_explore.push(other_end);
					}
				}
			}

			cur_component_id++;
		}
	}

	return components_node_ids;
}

pair<vector<ulong>,vector<set<ulong>>> PolyLineGraph2D::compute_components() {
	vector<ulong> components_node_ids(get_nodes_amount());
	vector<set<ulong>> components;

	ulong other_end;

	bool* explored = new bool[get_nodes_amount()];
	for(ulong i=0; i < get_nodes_amount(); i++)
		explored[i] = false;

	bool* in_to_explore = new bool[get_nodes_amount()];
	for(ulong i=0; i < get_nodes_amount(); i++)
		in_to_explore[i] = false;

	stack<ulong> to_explore;

	ulong cur_component_id = 0;

	// Explore
	for(ulong start_node = 0; start_node < get_nodes_amount(); start_node++) {
		if(!explored[start_node]) {
			vector<glm::vec4> cur_res;
			explored[start_node] = true;

			set<ulong> cur_component;
			cur_component.insert(start_node);
			components_node_ids[start_node] = cur_component_id;

			for(const auto p_id : connections[start_node]) {
				other_end = polylines[p_id].get_other_end(start_node);
				in_to_explore[other_end] = true;
				to_explore.push(other_end);
			}

			while(!to_explore.empty()) {
				const ulong cur_node = to_explore.top();
				to_explore.pop();

				cur_component.insert(cur_node);
				components_node_ids[cur_node] = cur_component_id;

				in_to_explore[cur_node] = false;
				explored[cur_node] = true;

				for(const auto p_id : connections[cur_node]) {
					other_end = polylines[p_id].get_other_end(cur_node);
					// " || other_end == cur_node " deals with 1-loops
					if(!explored[other_end] || other_end == cur_node) {
						if(!in_to_explore[other_end] && other_end != cur_node)
							to_explore.push(other_end);
					}
				}
			}

			components.push_back(cur_component);
			cur_component_id++;
		}
	}

	return make_pair(components_node_ids,components);
}

PolyLineGraph2D::plg_components PolyLineGraph2D::compute_plg_components()
{
	pair<pair<vector<ulong>,vector<set<ulong>>>,pair<vector<ulong>,vector<set<ulong>>>> tmp = compute_components_with_polylines();
	return plg_components(tmp.first.first,tmp.first.second, tmp.second.first,tmp.second.second);
}

pair<pair<vector<ulong>,vector<set<ulong>>>,pair<vector<ulong>,vector<set<ulong>>>> PolyLineGraph2D::compute_components_with_polylines() {
	pair<vector<ulong>,vector<set<ulong>>> component_nodes = compute_components();
	vector<ulong> polyline_component_ids(get_polylines_amount());
	vector<set<ulong>> polyline_components(component_nodes.second.size());
	for(ulong i=0; i < get_nodes_amount(); i++) {
		if(is_valid_node(i)) {
			ulong component_id = component_nodes.first[i];
			for(auto &pid : connections[i]) {
				polyline_component_ids[pid] = component_id;
				polyline_components[component_id].insert(pid);
			}
		}
	}

	return make_pair(component_nodes, make_pair(polyline_component_ids,polyline_components));
}

void PolyLineGraph2D::filter_components_by_polylinesmoothlength() {
	pair<pair<vector<ulong>,vector<set<ulong>>>,pair<vector<ulong>,vector<set<ulong>>>> components_data = compute_components_with_polylines();
	vector<float> polyline_smoothlengths(get_polylines_amount());
	for(ulong i=0; i < get_polylines_amount(); i++)
		if(is_valid_polyline(i))
			polyline_smoothlengths[i] = polylines[i].compute_max_smooth_length();


/*
	std::sort(polyline_smoothlengths.begin(),polyline_smoothlengths.end());
	// smoothlength_filter
	const float smoothlength_filter = polyline_smoothlengths[polyline_smoothlengths.size() * 0.8];
*/

	vector<float> polyline_smoothlengths_cpy = polyline_smoothlengths;
	const ulong smoothlength_filter_index = polyline_smoothlengths.size() * TOP_FILTER_BY_POLYLINESMOOTHLENGTH;
	std::nth_element (polyline_smoothlengths.begin(), polyline_smoothlengths.begin()+smoothlength_filter_index, polyline_smoothlengths.end());
	const float smoothlength_filter = polyline_smoothlengths[smoothlength_filter_index];

	// Filter components without any polyline above smoothlength_filter
	const vector<ulong> &polyline_components_ids = components_data.second.first;
	const vector<set<ulong>> &polyline_components = components_data.second.second;
	vector<bool> components_to_remove(polyline_components.size());
	for(ulong i=0; i < components_to_remove.size(); i++)
		components_to_remove[i] = true;


	for(ulong component_id = 0; component_id < polyline_components.size(); component_id++)
	{
		for(auto &pid : polyline_components[component_id])
			if(polyline_smoothlengths_cpy[pid] >= smoothlength_filter) {
				components_to_remove[component_id] = false;
				break;
			}
	}

	for(ulong i=0; i < components_to_remove.size(); i++)
		if(components_to_remove[i])
			for(set<ulong>::iterator it = polyline_components[i].begin(); it != polyline_components[i].end(); it++)
				remove_polyline(*it);

}

vector<pair<ulong,vector<PolyLineGraph2D::polyline::pl_point>>> PolyLineGraph2D::intersect_polylines(const glm::vec4 &segment)
{
	vector<pair<ulong,vector<PolyLineGraph2D::polyline::pl_point>>> res;

	for(ulong i=0; i < polylines.size();i++)
		if(is_valid_polyline(i)) {
			vector<PolyLineGraph2D::polyline::pl_point> cur_res = polylines[i].intersect_segment(segment);
			if(cur_res.size() > 0)
				res.push_back(make_pair(i,cur_res));
		}

	return res;
}

void PolyLineGraph2D::intersect_ray_first_polyline_within_dist(const ray2d &r, const float max_distance, plg_point &plgp, float &distance, bool &found) {
	polyline::pl_point closest_plp,cur_plp;
	ulong closest_pl_id;
	bool intersection_found;
	distance = max_distance+1;
	float cur_dist;

	for(ulong i=0; i < polylines.size();i++)
		if(is_valid_polyline(i)) {
			polylines[i].first_intersect_ray(r,cur_plp,cur_dist,intersection_found);
			if(intersection_found && cur_dist < distance) {
				found = true;
				closest_plp = cur_plp;
				distance = cur_dist;
				closest_pl_id = i;
			}
		}

	if(found)
		plgp = plg_point(closest_pl_id, closest_plp);
}


void PolyLineGraph2D::intersect_ray_first_polyline_within_dist(const ray2d &r, const ulong polyline_to_exclude, const float max_distance, plg_point &plgp, float &distance, bool &found) {
	polyline::pl_point closest_plp,cur_plp;
	ulong closest_pl_id;
	bool intersection_found;
	distance = max_distance;
	float cur_dist;

	found = false;

	for(ulong i=0; i < polylines.size();i++)
		if(is_valid_polyline(i) && i!= polyline_to_exclude) {
			//polylines[i].first_intersect_ray(r,cur_plp,cur_dist,intersection_found);
			polylines[i].first_intersect_ray_approx(r,cur_plp,cur_dist,intersection_found);
			if(intersection_found && cur_dist <= distance) {
				found = true;
				closest_plp = cur_plp;
				distance = cur_dist;
				closest_pl_id = i;
			}
		}

	if(found)
		plgp = plg_point(closest_pl_id, closest_plp);
}

float PolyLineGraph2D::cpf_find_unbound(const glm::vec2 &coords, PolyLineGraph2D::plg_point &plgp) {
	float min_dist = cpf_find_unbound(coords,plgp.polyline_id,plgp.plp.segment_index,plgp.plp.coords);
	return min_dist;
}

float PolyLineGraph2D::cpf_find_unbound(const glm::vec2 &coords, ulong &closest_polyline, ulong &closest_segment_on_polyline, glm::vec2 &projection_coords) {
	float min_distsq,cur_dist;
	ulong cur_segm;
	glm::vec2 cur_prj;

	min_distsq = std::numeric_limits<float>::max();

	for(ulong polyline_id = 0; polyline_id < polylines.size(); polyline_id++)
		if(is_valid_polyline(polyline_id)) {
			const polyline &p = polylines[polyline_id];
			cur_dist = p.compute_distancesq(coords,cur_segm,cur_prj);
			if(cur_dist < min_distsq) {
				min_distsq = cur_dist;
				closest_polyline = polyline_id;
				closest_segment_on_polyline = cur_segm;
				projection_coords = cur_prj;
			}
		}

	return sqrt(min_distsq);
}

vector<tuple<ulong,ulong,glm::vec2,float>> PolyLineGraph2D::cpf_find_within_radius(const glm::vec2 &coords, const float max_dist) {
	float cur_distsq;
	ulong cur_segm;
	glm::vec2 cur_prj;
	const float max_dist_sq = max_dist * max_dist;

	vector<tuple<ulong,ulong,glm::vec2,float>> res;

	for(ulong polyline_id = 0; polyline_id < polylines.size(); polyline_id++)
		if(is_valid_polyline(polyline_id)) {
			const polyline &p = polylines[polyline_id];
			cur_distsq = p.compute_distancesq(coords,cur_segm,cur_prj);
			if(cur_distsq < max_dist_sq) {
				res.push_back(make_tuple(polyline_id,cur_segm,cur_prj,sqrt(cur_distsq)));
			}
		}

	return res;
}

bool one_or_zero_correspondences(const pair<PolyLineGraph2D::plg_point,std::vector<std::vector<PolyLineGraph2D::plg_point>>> &p)
{

	for(const auto &v : p.second)
		if(v.size() > 1)
			return false;

	return true;
}

int amount_of_total_2d_correspondences(const pair<PolyLineGraph2D::plg_point,std::vector<std::vector<PolyLineGraph2D::plg_point>>> &p)
{
	int count = 0;
	for(const auto &v : p.second)
		count += v.size();
	return count;
}

vector<glm::vec2> convert_vecplgp_to_vec2(const std::vector<PolyLineGraph2D::plg_point> &plgps)
{
	vector<glm::vec2> coords;
	for(const auto &plgp : plgps)
		coords.push_back(plgp.plp.coords);
	return coords;
}

pair<vector<glm::vec2>,vector<int>> convert_plgpoint_correspondences(const vector<int> &cams, const std::vector<std::vector<PolyLineGraph2D::plg_point>> &potential_match_correspondences)
{
	vector<glm::vec2> coords;
	vector<int> ids;

	for(int i=0; i < potential_match_correspondences.size(); i++) {
		const int cam_id = cams[i];
		const std::vector<PolyLineGraph2D::plg_point> &img_potential_match = potential_match_correspondences[i];
		if(img_potential_match.size() > 0)
		{
			ids.push_back(cam_id);
			coords.push_back(img_potential_match[0].plp.coords);
		}
	}

	return make_pair(coords,ids);
}

pair<vector<PolyLineGraph2D::plg_point>,vector<int>> convert_plgpoint_correspondences_plgp(const vector<int> &cams, const std::vector<std::vector<PolyLineGraph2D::plg_point>> &potential_match_correspondences)
{
	vector<PolyLineGraph2D::plg_point> plgps;
	vector<int> ids;

	for(int i=0; i < potential_match_correspondences.size(); i++) {
		const int cam_id = cams[i];
		const std::vector<PolyLineGraph2D::plg_point> &img_potential_match = potential_match_correspondences[i];
		if(img_potential_match.size() > 0)
		{
			ids.push_back(cam_id);
			plgps.push_back(img_potential_match[0]);
		}
	}

	return make_pair(plgps,ids);
}
