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


#ifndef INCLUDE_EDGEGRAPH3D_PLGS_GRAPH_EDGE_ADJACENCY_LIST_HPP_
#define INCLUDE_EDGEGRAPH3D_PLGS_GRAPH_EDGE_ADJACENCY_LIST_HPP_

#include <pair>
#include <utility>
#include <vector>

#include "graph.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"

template <typename T>
class GraphAdjacencyEdgeList : public Graph<T> {
	/**
	 * hub nodes data:
	 * 		vector<multi_edge>
	 * intermediate nodes:
	 * 		not existing in graph structure
	 *
	 */
public:
	/**
	 * start is not included.
	 *
	 */
	struct edge {
		std::pair<ulong> extreme_nodes;
		vector<ulong> edge_nodes;
	};

private:
	vector<vector<ulong>> adjacency_lists;
protected:
	void increase_size(ulong sz);
public:
	GraphAdjacencyEdgeList();
	GraphAdjacencyEdgeList(ulong input_nodes_num);

	void add_edge(ulong start_node,ulong end_node,vector<ulong> intermediate_nodes);

	vector<vector<ulong>> get_adjacency_lists();

};

template <typename T>
GraphAdjacencyEdgeList<T>::GraphAdjacencyEdgeList() : Graph<T>() {}

template <typename T>
GraphAdjacencyEdgeList<T>::GraphAdjacencyEdgeList(ulong input_nodes_num) : Graph<T>(input_nodes_num) {
	for(int i=0; i<input_nodes_num; i++)
		adjacency_lists.push_back(vector<ulong>());
}

template <typename T>
void GraphAdjacencyEdgeList<T>::increase_size(ulong sz) {
	if(sz > get_nodes_num()) {
		adjacency_lists.resize(sz);
		for(int i=get_nodes_num(); i < sz; i++) {
			adjacency_lists[i] = vector<ulong>();
		}
	}
	Graph<T>::increase_size(sz);
}

template <typename T>
void GraphAdjacencyEdgeList<T>::add_edge(ulong start_node,ulong end_node) {
	if(start_node >= get_nodes_num())
		increase_size(start_node+1);
	if(end_node >= get_nodes_num())
		increase_size(end_node+1);
	adjacency_lists[start_node].push_back(end_node);
}

template <typename T>
vector<vector<ulong>> GraphAdjacencyEdgeList<T>::get_adjacency_lists() {
	return adjacency_lists;
}



#endif /* INCLUDE_EDGE_MATCHER_EDGES_OUTPUT_CONVERTER_GRAPHADJACENCYLIST_HPP_ */
