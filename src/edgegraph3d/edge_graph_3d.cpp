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

#include <iostream>

#include <omp.h>

#include "edge_matcher.hpp"
#include "edge_graph_3d_input_params.hpp"

int main(int argc, char** argv) {
	edge_matcher_input_params emip;

	if(!read_edge_matcher_input_params(argc,argv,emip))
		return 1;

	std::cout << "\n\nEdgeGraph3D\n===========\n\n";

	std::cout << "Input SfM data: " << emip.sfm_data_file << "\n";
	std::cout << "Input images: " << emip.images_folder << "\n";
	std::cout << "Input edges: " << emip.input_edges_folder << "\n";
	std::cout << "Working folder: " << emip.em_out_folder << "\n";
	std::cout << "Output: " << emip.output_json << "\n";

	std::cout << "\n";

#ifdef SWITCH_RUNPARALLEL
	std::cout << "Computing on " << omp_get_max_threads() << " threads\n\n";
#endif

	edge_matching(emip);

	std::cout << "\nDone.\n\n";
}
