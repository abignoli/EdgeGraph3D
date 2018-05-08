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

#include "community_detection_interface.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "defs.h"
#include "driverForGraphClustering_edited.hpp"
#include "graph_adjacency_set_undirected_no_type_weighted.hpp"

// Read Grappolo's output communities
vector<long> read_cluster_info(const char *input_filename) {
	vector<long> res;
	  string cur_line;
	  ifstream input_file (input_filename);
	  if (input_file.is_open())
	  {
	    while ( getline (input_file,cur_line) )
	    	res.push_back(stoul(cur_line,nullptr,10));
	    input_file.close();
	  } else
		  cout << "Unable to open file";

	  return res;
}

vector<long> compute_communities(const GraphAdjacencySetUndirectedNoTypeWeighted &pmgw, const char *out_compatibility_graph_file, const char *out_polyline_communities) {
	pmgw.write_to_file(out_compatibility_graph_file);

	// Create options for Grappolo's community detection
	clustering_parameters opts;
	opts.inFile = out_compatibility_graph_file;
	opts.ftype = 2;
	opts.coloring = true;
	opts.output = true;
	opts.minGraphSize = 10;

	// Compute communities
	grappolo_community_detection(opts, out_polyline_communities);

	// Read output communities
	return read_cluster_info(out_polyline_communities);
}


