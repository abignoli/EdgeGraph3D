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


#ifndef INCLUDE_EDGEGRAPH3D_PLGS_GRAPH_HPP_
#define INCLUDE_EDGEGRAPH3D_PLGS_GRAPH_HPP_

#include <vector>

#include "edge_graph_3d_utilities.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"

using namespace std;

template <typename T>
class Graph {
private:
	vector<T> nodes_data;
	ulong nodes_num;
protected:
	void increase_size(ulong sz);
public:
	Graph();
	Graph(ulong input_nodes_num);
	T get_node_data(ulong id);
	ulong get_nodes_num();
	void set_node_data(ulong id, T data);
	vector<T> get_nodes_data();
	virtual void add_edge(ulong start_node,ulong end_node)=0;
};

template <typename T>
Graph<T>::Graph() {
	nodes_num = 0;
}

template <typename T>
void Graph<T>::increase_size(ulong sz) {
	if(sz > nodes_num) {
		nodes_data.resize(sz);
		nodes_num=sz;
	}
}

template <typename T>
ulong Graph<T>::get_nodes_num() {
	return nodes_num;
}

template <typename T>
void Graph<T>::set_node_data(ulong id, T data) {
	nodes_data[id] = data;
}

template <typename T>
T Graph<T>::get_node_data(ulong id) {
	return nodes_data[id];
}

template <typename T>
Graph<T>::Graph(ulong input_nodes_num) : nodes_num(input_nodes_num) {
	nodes_data.resize(input_nodes_num);
}

template <typename T>
vector<T> Graph<T>::get_nodes_data() {
	return nodes_data;
}

#endif /* INCLUDE_EDGEGRAPH3D_PLGS_GRAPH_HPP_ */
