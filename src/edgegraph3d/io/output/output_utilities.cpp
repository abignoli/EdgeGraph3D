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

#include "output_utilities.hpp"

#include <opencv2/core/mat.hpp>
#include <opencv2/core/mat.inl.hpp>
#include <fstream>
#include <string>
#include <vector>

#include "polyline_graph_2d_hmap_impl.hpp"
#include "SfMData.h"
#include "filtering_close_plgps.hpp"
#include "output_sfm_data.hpp"
#include "polyline_graph_2d.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "glm.hpp"


#define FEATURE_ID 0
void print_new_point_to_file(const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &p3d, ofstream &outfile, const ulong current_point_id, const ulong last_point_id) {
	outfile <<
"        {" << endl <<
"            \"key\": " << current_point_id << "," << endl <<
"            \"value\": {" << endl <<
"                \"X\": [" << endl <<
"                    " << get<0>(p3d).x << "," << endl <<
"                    " << get<0>(p3d).y << "," << endl <<
"                    " << get<0>(p3d).z << endl <<
"                ],\n" << endl <<
"                \"observations\": [" << endl;
	for(int i=0; i < get<1>(p3d).size(); i++) {
		outfile <<
"                    {" << endl <<
"                        \"key\": " << get<2>(p3d)[i] << "," << endl <<
"                        \"value\": {" << endl <<
"                            \"id_feat\": " << FEATURE_ID << "," << endl <<
"                            \"x\": [" << endl <<
"                                " << get<1>(p3d)[i].plp.coords.x << "," << endl <<
"                                " << get<1>(p3d)[i].plp.coords.y << endl <<
"                            ]" << endl <<
"                        }" << endl <<
"                    }";
		if(i < get<1>(p3d).size()-1)
			outfile << ",";
		outfile << endl;
	}
	outfile <<
"                ]" << endl <<
"            }" << endl <<
"        }";
	if(current_point_id!=last_point_id)
		outfile << ",";
	outfile << endl;
}

void print_new_points_to_file(const SfMData &sfmd, const vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &p3ds, const string &outfilename) {
	ofstream outfile;
	outfile.open(outfilename);
	ulong current_point_id = sfmd.numPoints_;
	ulong last_point_id = sfmd.numPoints_ + p3ds.size() - 1;
	for(const auto &p3d:p3ds) {
		print_new_point_to_file(p3d,outfile,current_point_id,last_point_id);
		current_point_id++;
	}
	outfile.close();
}

void add_3dpoints_to_sfmd(SfMData &sfmd, const vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &p3ds) {
	//print_new_points_to_file(sfmd,p3ds,newpoints_outfilename);
	for(const auto &p3d: p3ds) {
		int new_point_id = sfmd.points_.size();

		sfmd.points_.push_back(get<0>(p3d));
		for(const auto cam_id : get<2>(p3d))
			sfmd.pointsVisibleFromCamN_[cam_id].push_back(new_point_id);

		vector<glm::vec2> coords = convert_vecplgp_to_vec2(get<1>(p3d));
		sfmd.point2DoncamViewingPoint_.push_back(coords);

		sfmd.camViewingPointN_.push_back(vector<int>(get<2>(p3d).begin(), get<2>(p3d).end()));
	}
	sfmd.numPoints_ = sfmd.points_.size();
}

void write_sfmd_plus_points(const char *original_sfm_data_file, vector<Mat> &imgs, const SfMData &sfmd, const vector<PolyLineGraph2DHMapImpl> &plgs, const vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &p3ds, const char *out_folder, const char *newSfMname) {
	SfMData sfmd_new(sfmd);
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> p3ds_new = filter_3d_points_close_2d_array(plgs,imgs[0].size(),p3ds);
	add_3dpoints_to_sfmd(sfmd_new,p3ds_new);
	output_sfm_data(original_sfm_data_file, sfmd_new, out_folder + std::string(newSfMname));

}

