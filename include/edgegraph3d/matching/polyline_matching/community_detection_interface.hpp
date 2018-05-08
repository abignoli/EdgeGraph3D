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

#ifndef INCLUDE_EDGEGRAPH3D_MATCHING_POLYLINE_MATCHING_COMMUNITY_DETECTION_INTERFACE_HPP_
#define INCLUDE_EDGEGRAPH3D_MATCHING_POLYLINE_MATCHING_COMMUNITY_DETECTION_INTERFACE_HPP_

#include <vector>

#include "graph_adjacency_set_undirected_no_type_weighted.hpp"

vector<long> compute_communities(const GraphAdjacencySetUndirectedNoTypeWeighted &pmgw, const char *out_compatibility_graph_file, const char *out_polyline_communities);

#endif /* INCLUDE_EDGEGRAPH3D_MATCHING_POLYLINE_MATCHING_COMMUNITY_DETECTION_INTERFACE_HPP_ */
