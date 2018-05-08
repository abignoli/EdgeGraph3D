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


#ifndef INCLUDE_EDGEGRAPH3D_PLGS_GRAPH_ADJACENCY_SET_NO_TYPE_HPP_
#define INCLUDE_EDGEGRAPH3D_PLGS_GRAPH_ADJACENCY_SET_NO_TYPE_HPP_

#include <set>
#include <vector>

#include "graph_no_type.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"

class GraphAdjacencySetNoType : public GraphNoType {
protected:
	vector<bool> visited; // Just used for internal computations
	void increase_size(ulong sz);
public:
	vector<std::set<ulong>> adjacency_lists;
	ulong get_edges_amount();
	ulong add_node();
	GraphAdjacencySetNoType();
	GraphAdjacencySetNoType(ulong input_nodes_num);
	void add_edge(ulong start_node,ulong end_node);
	vector<std::set<ulong>> get_adjacency_sets();
	bool is_connected(ulong start_node,ulong end_node);
	bool is_connected(ulong start_node,ulong end_node,ulong max_dist);
};




#endif /* INCLUDE_EDGE_MATCHER_EDGES_OUTPUT_CONVERTER_GRAPHADJACENCYLIST_HPP_ */
