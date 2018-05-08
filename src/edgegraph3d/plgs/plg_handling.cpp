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


#include "plg_handling.hpp"

#include <opencv2/core/cvstd.hpp>
#include <opencv2/core/mat.hpp>
#include <opencv2/core/matx.hpp>
#include <iostream>
#include <string>
#include <vector>

#include "polyline_graph_2d_hmap_impl.hpp"
#include "SfMData.h"
#include "convert_edge_images_pixel_to_segment.hpp"
#include "edge_graph_3d_utilities.hpp"

void write_plgs(const vector<PolyLineGraph2DHMapImpl> &plgs, const SfMData &sfmd, char *em_out_folder,  const String &s) {
	for(int i=0; i < plgs.size(); i++) {
		const PolyLineGraph2DHMapImpl &plg = plgs[i];
		cout << "Writing " << compute_plg_path(em_out_folder, sfmd, i,s) << "\n";
		serialize_plg(plg,compute_plg_path(em_out_folder, sfmd, i,s));
	}
}

vector<PolyLineGraph2DHMapImpl> compute_and_write_plgs(const vector<Mat> &edge_imgs, const Vec<unsigned char, 3> &edge_color, const SfMData &sfmd, char *em_out_folder,  const String &s) {
	vector<PolyLineGraph2DHMapImpl> plgs = convert_edge_images_to_optimized_polyline_graphs(edge_imgs,edge_color);
	write_plgs(plgs, sfmd, em_out_folder, s);
	return plgs;
}

vector<PolyLineGraph2DHMapImpl> read_plgs(const SfMData &sfmd, char *em_out_folder,  const String &s) {
	vector<PolyLineGraph2DHMapImpl> res;

	for(int i=0; i < sfmd.numCameras_; i++) {
		res.push_back(deserialize_plg(compute_plg_path(em_out_folder, sfmd, i,s)));
	}

	return res;
}

// This method is mostly used for experimentation and is not referenced in production code
void transform_plgs(vector<PolyLineGraph2DHMapImpl> &plgs)
{
	for(auto &plg : plgs) {
//		plg.filter_components_by_polylinesmoothlength();
//		plg.connect_close_extremes();
//		plg.remove_2connection_nodes();
//		plg.simplify();

		cout << "Optimizing PLG...\n";

		plg.optimize();

//		plg.filter_components_by_polylinesmoothlength();
//		plg.connect_close_extremes_following_direction();
//		plg.connect_close_extremes();
//		plg.remove_2connection_nodes();
//		plg.simplify();
		plg.connect_close_extremes();
		plg.remove_2connection_nodes();
		plg.simplify();

	}

}



