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

#include <fstream>
#include <iostream>
#include "graph_adjacency_set_undirected_no_type_weighted.hpp"

GraphAdjacencySetUndirectedNoTypeWeighted::GraphAdjacencySetUndirectedNoTypeWeighted() : GraphAdjacencySetUndirectedNoType() {}

GraphAdjacencySetUndirectedNoTypeWeighted::GraphAdjacencySetUndirectedNoTypeWeighted(ulong input_nodes_num) : GraphAdjacencySetUndirectedNoType(input_nodes_num) {}

void GraphAdjacencySetUndirectedNoTypeWeighted::add_edge(ulong node_a,ulong node_b, float weight) {
	weights[make_pair(node_a,node_b)] = weight;
	GraphAdjacencySetUndirectedNoType::add_edge(node_a,node_b);
}

float GraphAdjacencySetUndirectedNoTypeWeighted::get_weight(ulong node_a,ulong node_b, bool &found) {
	undirectedweightmap::const_iterator it = weights.find(make_pair(node_a,node_b));
	found = false;

	if(it != weights.end()) {
		found = true;
		return it->second;
	}
	return 0.0;
}

void GraphAdjacencySetUndirectedNoTypeWeighted::write_to_file(char *outputfile) {
	  std::ofstream file;
	  file.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
	  try {
	    file.open (outputfile);

	    ulong amount_of_nodes = get_nodes_num();
	    ulong amount_of_edges = get_edges_amount();

	    file << "p sp " << amount_of_nodes << " " << amount_of_edges << "\n";
		for(ulong node1=0; node1 < get_nodes_num(); node1++)
			for(const auto node2: adjacency_lists[node1])
				file << "a " << node1+1 << " " << node2+1 << " " << ((weights.find(make_pair(node1,node2)))->second) << "\n";

	    file.close();
	  }
	  catch (std::ifstream::failure e) {
	    std::cerr << "Exception opening graph file for writing \n";
	  }

}
