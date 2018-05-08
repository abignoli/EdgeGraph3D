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


#ifndef INCLUDE_EDGEGRAPH3D_PLGS_GRAPH_ADJACENCY_LIST_HPP_
#define INCLUDE_EDGEGRAPH3D_PLGS_GRAPH_ADJACENCY_LIST_HPP_

#include <vector>

#include "graph.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"

template <typename T>
class GraphAdjacencyList : public Graph<T> {
private:
	vector<vector<ulong>> adjacency_lists;
protected:
	void increase_size(ulong sz);
public:
	GraphAdjacencyList();
	GraphAdjacencyList(ulong input_nodes_num);
	void add_edge(ulong start_node,ulong end_node);
	vector<vector<ulong>> get_adjacency_lists();
};

template <typename T>
GraphAdjacencyList<T>::GraphAdjacencyList() : Graph<T>() {}

template <typename T>
GraphAdjacencyList<T>::GraphAdjacencyList(ulong input_nodes_num) : Graph<T>(input_nodes_num) {
	for(int i=0; i<input_nodes_num; i++)
		adjacency_lists.push_back(vector<ulong>());
}

template <typename T>
void GraphAdjacencyList<T>::increase_size(ulong sz) {
	if(sz > get_nodes_num()) {
		adjacency_lists.resize(sz);
		for(int i=get_nodes_num(); i < sz; i++) {
			adjacency_lists[i] = vector<ulong>();
		}
	}
	Graph<T>::increase_size(sz);
}

template <typename T>
void GraphAdjacencyList<T>::add_edge(ulong start_node,ulong end_node) {
	if(start_node >= get_nodes_num())
		increase_size(start_node+1);
	if(end_node >= get_nodes_num())
		increase_size(end_node+1);
	adjacency_lists[start_node].push_back(end_node);
}

template <typename T>
vector<vector<ulong>> GraphAdjacencyList<T>::get_adjacency_lists() {
	return adjacency_lists;
}



#endif /* INCLUDE_EDGEGRAPH3D_PLGS_GRAPH_ADJACENCY_LIST_HPP_ */
