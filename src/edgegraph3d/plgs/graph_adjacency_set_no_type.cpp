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

#include "graph_adjacency_set_no_type.hpp"

#include <map>
#include <stack>

GraphAdjacencySetNoType::GraphAdjacencySetNoType() : GraphNoType() {}


GraphAdjacencySetNoType::GraphAdjacencySetNoType(ulong input_nodes_num) : GraphNoType(input_nodes_num) {
	for(int i=0; i<input_nodes_num; i++) {
		adjacency_lists.push_back(std::set<ulong>());
		visited.push_back(false);
	}
}


void GraphAdjacencySetNoType::increase_size(ulong sz) {
	if(sz > get_nodes_num()) {
		ulong cur_sz = adjacency_lists.size();
		adjacency_lists.resize(sz);
		for(ulong i=cur_sz; i < sz; i++) {
			adjacency_lists[i] = std::set<ulong>();
			visited.push_back(false);
		}
	}
	GraphNoType::increase_size(sz);
}

ulong GraphAdjacencySetNoType::get_edges_amount() {
	ulong count = 0;

	for(const auto &s : adjacency_lists)
		count += s.size();

	return count;
}

ulong GraphAdjacencySetNoType::add_node() {
	adjacency_lists.push_back(std::set<ulong>());
	visited.push_back(false);
	return GraphNoType::add_node();
}


void GraphAdjacencySetNoType::add_edge(ulong start_node,ulong end_node) {
	if(start_node >= get_nodes_num())
		increase_size(start_node+1);
	if(end_node >= get_nodes_num())
		increase_size(end_node+1);
	adjacency_lists[start_node].insert(end_node);
}


vector<std::set<ulong>> GraphAdjacencySetNoType::get_adjacency_sets() {
	return adjacency_lists;
}

bool connections_contains_node(const set<ulong> &node_connections, ulong other_end) {
	return node_connections.find(other_end) != node_connections.end();
}

bool GraphAdjacencySetNoType::is_connected(ulong start_node,ulong end_node) {
	return connections_contains_node(adjacency_lists[start_node], end_node);
}

bool GraphAdjacencySetNoType::is_connected(ulong start_node,ulong end_node, ulong max_dist) {
	vector<bool> visited_vec;
	visited_vec.push_back(start_node);
	visited[start_node]=true;

/*	set<ulong> visited;
	visited.insert(start_node);*/
	bool found = false;

	stack<ulong> cur_to_visit;
	cur_to_visit.push(start_node);
	ulong cur_dist = 0;

	while(cur_dist <= max_dist && !cur_to_visit.empty()) {
		stack<ulong> next_to_visit;

		while(!found && !cur_to_visit.empty()) {
			const ulong cur_node = cur_to_visit.top();
			cur_to_visit.pop();

			for(const auto connected_node : adjacency_lists[cur_node]) {
				if(connected_node == end_node) {
					found = true;
					break;
				}

				if(!visited[connected_node]) {
					next_to_visit.push(connected_node);
					visited[connected_node] = true;
					visited_vec.push_back(connected_node);
				}
			}
		}

		for(auto v : visited_vec)
			visited[v] = false;

		cur_to_visit = next_to_visit;
		cur_dist++;
	}

	return found;
}
