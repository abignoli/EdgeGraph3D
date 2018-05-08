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


#ifndef INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_POLYLINE_2D_MAP_SEARCH_HPP_
#define INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_POLYLINE_2D_MAP_SEARCH_HPP_

#include <opencv2/core/types.hpp>
#include <set>
#include <vector>

#include "polyline_graph_2d.hpp"
#include "polyline_2d_map.hpp"
#include "polyline_graph_2d_hmap_impl.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"
#include "global_switches.hpp"
#include "glm.hpp"

class PolyLine2DMapSearch : PolyLine2DMap {
public:
	const float search_dist;
	const float search_dist_sq;

	set<ulong> PolyLine2DMapSearch::find_polylines_potentially_within_search_dist(const glm::vec2 &coords);
	void find_unique_polyline_potentially_within_search_dist(const glm::vec2 &coords, ulong &pl_id, bool &valid);
	vector<ulong> find_polylines_within_search_dist(const glm::vec2 &coords);
	vector<ulong> find_polylines_within_smaller_search_dist(const glm::vec2 &coords, const float smaller_search_dist);
	vector<tuple<ulong,ulong,glm::vec2,float>> find_polylines_within_search_dist_with_reprojections(const glm::vec2 &coords);
	std::vector<PolyLineGraph2D::plg_point> find_polylines_within_search_dist_with_plgps(const glm::vec2 &coords);
	std::vector<PolyLineGraph2D::plg_point> find_polylines_within_smaller_search_dist_with_plgps(const glm::vec2 &coords, const float smaller_search_dist);

	PolyLine2DMapSearch(const PolyLineGraph2DHMapImpl &plg, const Size img_sz, const float search_dist);

	~PolyLine2DMapSearch();
};


#endif /* INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_POLYLINE_2D_MAP_SEARCH_HPP_ */
