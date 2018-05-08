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


#ifndef INCLUDE_EDGEGRAPH3D_PLGS_GRAPH_NO_TYPE_HPP_
#define INCLUDE_EDGEGRAPH3D_PLGS_GRAPH_NO_TYPE_HPP_

#include "edge_graph_3d_utilities.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"

using namespace std;

class GraphNoType {
private:
	ulong nodes_num;
protected:
	void increase_size(ulong sz);
public:
	ulong add_node(); // Returns new node ID
	GraphNoType();
	GraphNoType(ulong input_nodes_num);
	ulong get_nodes_num();
	virtual void add_edge(ulong start_node,ulong end_node)=0;
};



#endif /* INCLUDE_EDGE_MATCHER_EDGES_OUTPUT_CONVERTER_GRAPH_HPP_ */
