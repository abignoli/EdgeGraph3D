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


#ifndef INCLUDE_EDGEGRAPH3D_PLGS_GRAPH_ADJACENCY_LIST_UNDIRECTED_HPP_
#define INCLUDE_EDGEGRAPH3D_PLGS_GRAPH_ADJACENCY_LIST_UNDIRECTED_HPP_

#include "graph_adjacency_list.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"

template <typename T>
class GraphAdjacencyListUndirected : public GraphAdjacencyList<T> {
public:
	GraphAdjacencyListUndirected();
	GraphAdjacencyListUndirected(ulong input_nodes_num);
	void add_edge(ulong node_a,ulong node_b);
};

template <typename T>
GraphAdjacencyListUndirected<T>::GraphAdjacencyListUndirected() : GraphAdjacencyList<T>() { }

template <typename T>
GraphAdjacencyListUndirected<T>::GraphAdjacencyListUndirected(ulong input_nodes_num) : GraphAdjacencyList<T>(input_nodes_num) { }

template <typename T>
void GraphAdjacencyListUndirected<T>::add_edge(ulong start_node,ulong end_node) {
	GraphAdjacencyList<T>::add_edge(start_node,end_node);
	GraphAdjacencyList<T>::add_edge(end_node,start_node);
}



#endif /* INCLUDE_EDGEGRAPH3D_PLGS_GRAPH_ADJACENCY_LIST_UNDIRECTED_HPP_ */
