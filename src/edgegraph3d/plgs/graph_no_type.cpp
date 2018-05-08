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

#include "graph_no_type.hpp"

GraphNoType::GraphNoType() {
	nodes_num = 0;
}

void GraphNoType::increase_size(ulong sz) {
	if(sz > nodes_num)
		nodes_num=sz;
}

ulong GraphNoType::add_node() {
	nodes_num++;
	return nodes_num-1;
}


ulong GraphNoType::get_nodes_num() {
	return nodes_num;
}



GraphNoType::GraphNoType(ulong input_nodes_num) : nodes_num(input_nodes_num) {}
