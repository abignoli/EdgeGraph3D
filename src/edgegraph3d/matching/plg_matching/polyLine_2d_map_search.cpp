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


#include "polyline_2d_map_search.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <map>
#include <stdexcept>
#include <utility>

#include "polyline_graph_2d_hmap_impl.hpp"
#include "glm.hpp"

PolyLine2DMapSearch::PolyLine2DMapSearch(const PolyLineGraph2DHMapImpl &plg, const Size img_sz, const float search_dist) :
	search_dist(search_dist), search_dist_sq(search_dist*search_dist), PolyLine2DMap(plg,img_sz,search_dist) {}

set<ulong> PolyLine2DMapSearch::find_polylines_potentially_within_search_dist(const glm::vec2 &coords) {
	set<ulong> res;
	bool on_boundary_row,on_boundary_col;

	if(coords.x <= 0 || coords.x >= img_sz.width || coords.y <= 0 || coords.y >= img_sz.height)
		return res;

	pair<ulong,ulong> cell_coords = get_2dmap_cell_from_coords(cell_dim, coords, on_boundary_row,on_boundary_col);

	cell_coords.first = cell_coords.first >= mapsz.width ? mapsz.width - 1 : cell_coords.first;
	cell_coords.second = cell_coords.second >= mapsz.height ? mapsz.height - 1 : cell_coords.second;

	if(!on_boundary_row && !on_boundary_col)
		// search 3x3 centered around cell coords
		for(int i= cell_coords.second > 0 ? -1 : 0; i <= (cell_coords.second < mapsz.height-1 ? 1 : 0); i++)
			for(int j= cell_coords.first > 0 ? -1 : 0; j<= (cell_coords.first < mapsz.width-1 ? 1 : 0); j++)
				std::copy( pls_id_maps[cell_coords.second+i][cell_coords.first+j].begin(), pls_id_maps[cell_coords.second+i][cell_coords.first+j].end(), std::inserter( res, res.end() ) );
	else if (on_boundary_row && !on_boundary_col)
		for(int i= cell_coords.second > 0 ? -1 : 0; i <= 0; i++)
			for(int j= cell_coords.first > 0 ? -1 : 0; j<= (cell_coords.first < mapsz.width-1 ? 1 : 0); j++)
				std::copy( pls_id_maps[cell_coords.second+i][cell_coords.first+j].begin(), pls_id_maps[cell_coords.second+i][cell_coords.first+j].end(), std::inserter( res, res.end() ) );
	else if (!on_boundary_row && on_boundary_col)
		for(int i= cell_coords.second > 0 ? -1 : 0; i <= (cell_coords.second < mapsz.height-1 ? 1 : 0); i++)
			for(int j= cell_coords.first > 0 ? -1 : 0; j<= 0; j++)
				std::copy( pls_id_maps[cell_coords.second+i][cell_coords.first+j].begin(), pls_id_maps[cell_coords.second+i][cell_coords.first+j].end(), std::inserter( res, res.end() ) );
	else if (on_boundary_row && on_boundary_col)
		for(int i= cell_coords.second > 0 ? -1 : 0; i <= 0; i++)
			for(int j= cell_coords.first > 0 ? -1 : 0; j<= 0; j++)
				std::copy( pls_id_maps[cell_coords.second+i][cell_coords.first+j].begin(), pls_id_maps[cell_coords.second+i][cell_coords.first+j].end(), std::inserter( res, res.end() ) );

	return res;
}

PolyLine2DMapSearch::~PolyLine2DMapSearch() {}

void PolyLine2DMapSearch::find_unique_polyline_potentially_within_search_dist(const glm::vec2 &coords, ulong &pl_id, bool &valid) {
	valid=false;
	set<ulong> tmp_res = find_polylines_potentially_within_search_dist(coords);
	if(tmp_res.size()==1) {
		valid = true;
		pl_id = *(tmp_res.begin());
	}
}

vector<ulong> PolyLine2DMapSearch::find_polylines_within_search_dist(const glm::vec2 &coords) {
	set<ulong> tmp_res = find_polylines_potentially_within_search_dist(coords);
	vector<ulong> res;

	ulong closest_segm;
	glm::vec2 projection;

	for(const auto pl_id : tmp_res)
		if(plg.polylines[pl_id].compute_distancesq(coords, closest_segm, projection) <= search_dist_sq)
			res.push_back(pl_id);

	return res;
}

vector<ulong> PolyLine2DMapSearch::find_polylines_within_smaller_search_dist(const glm::vec2 &coords, const float smaller_search_dist) {
	if(smaller_search_dist > search_dist)
		std::invalid_argument("find_polylines_within_smaller_search_dist - smaller_search_dist > search_dist");

	set<ulong> tmp_res = find_polylines_potentially_within_search_dist(coords);
	vector<ulong> res;

	ulong closest_segm;
	glm::vec2 projection;
	const float smaller_search_distsq = smaller_search_dist*smaller_search_dist;

	for(const auto pl_id : tmp_res)
		if(plg.polylines[pl_id].compute_distancesq(coords, closest_segm, projection) <= smaller_search_distsq)
			res.push_back(pl_id);

	return res;
}

vector<tuple<ulong,ulong,glm::vec2,float>> PolyLine2DMapSearch::find_polylines_within_search_dist_with_reprojections(const glm::vec2 &coords) {
	set<ulong> tmp_res = find_polylines_potentially_within_search_dist(coords);
	vector<tuple<ulong,ulong,glm::vec2,float>> res;

	ulong closest_segm;
	glm::vec2 projection;
	float cur_distsq;

	for(const auto pl_id : tmp_res) {
		cur_distsq = plg.polylines[pl_id].compute_distancesq(coords, closest_segm, projection);
		if(cur_distsq <= search_dist_sq)
			res.push_back(make_tuple(pl_id,closest_segm,projection,sqrt(cur_distsq)));
	}

	return res;
}

std::vector<PolyLineGraph2D::plg_point> PolyLine2DMapSearch::find_polylines_within_search_dist_with_plgps(const glm::vec2 &coords) {
	set<ulong> tmp_res = find_polylines_potentially_within_search_dist(coords);
	std::vector<PolyLineGraph2D::plg_point> res;

	ulong closest_segm;
	glm::vec2 projection;

	for(const auto pl_id : tmp_res)
		if(plg.polylines[pl_id].compute_distancesq(coords, closest_segm, projection) <= search_dist_sq)
			res.push_back(PolyLineGraph2D::plg_point(pl_id, closest_segm, projection));

	return res;
}

std::vector<PolyLineGraph2D::plg_point> PolyLine2DMapSearch::find_polylines_within_smaller_search_dist_with_plgps(const glm::vec2 &coords, const float smaller_search_dist) {
	if(smaller_search_dist > search_dist)
		std::invalid_argument("find_polylines_within_smaller_search_dist - smaller_search_dist > search_dist");

	set<ulong> tmp_res = find_polylines_potentially_within_search_dist(coords);
	std::vector<PolyLineGraph2D::plg_point> res;

	ulong closest_segm;
	glm::vec2 projection;

	const float smaller_search_distsq = smaller_search_dist*smaller_search_dist;

	for(const auto pl_id : tmp_res)
		if(plg.polylines[pl_id].compute_distancesq(coords, closest_segm, projection) <= smaller_search_distsq)
			res.push_back(PolyLineGraph2D::plg_point(pl_id, closest_segm, projection));

	return res;
}
