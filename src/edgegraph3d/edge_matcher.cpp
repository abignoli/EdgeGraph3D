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

#include "edge_matcher.hpp"

#include <omp.h>
#include <opencv2/core/mat.hpp>
#include <opencv2/core/mat.inl.hpp>
#include <opencv2/core/types.hpp>
#include <iostream>
#include <string>
#include <vector>

#include "polyline_graph_2d_hmap_impl.hpp"
#include "SfMData.h"
#include "plg_edge_manager.hpp"
#include "outliers_filtering.hpp"
#include "convert_edge_images_pixel_to_segment.hpp"
#include "edge_graph_3d_input_params.hpp"
#include "output_sfm_data.hpp"
#include "plgpcm_3views_plg_following.hpp"
#include "pipelines.hpp"
#include "plg_matches_manager.hpp"
#include "polyline_2d_map_search.hpp"
#include "polyline_graph_3d_hmap_impl.hpp"
#include "data_bundle.hpp"
#include "drawing_utilities.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "geometric_utilities.hpp"
#include "global_defines.hpp"

using namespace std;
using namespace cv;

int edge_matching(edge_matcher_input_params &emip) {
	SfMData sfmd = read_sfm_data(emip.sfm_data_file);
	return edge_matching(emip,sfmd);
}

int edge_matching(edge_matcher_input_params &emip, SfMData &sfm_data) {

	vector<Mat> imgs;

	if(parse_images(emip.images_folder, sfm_data, imgs))
		return -1; // something went wrong in loading

	vector<Mat> input_edge_imgs;

	if(parse_images(emip.input_edges_folder, sfm_data, input_edge_imgs))
		return -1; // something went wrong in loading

	Size img_size = Size(imgs[0].cols,imgs[0].rows);

	double plgs_time_start = omp_get_wtime();

	// Compute PolyLine Graphs from input edge_images
	vector<PolyLineGraph2DHMapImpl> plgs = convert_edge_images_to_optimized_polyline_graphs(input_edge_imgs,EDGE_COLOR); // no writing plgs

	double plgs_time_end = omp_get_wtime();
	double plgs_time = plgs_time_end - plgs_time_start;
	cout << "PLGS computation time: " << plgs_time << "s" << endl;

	if(emip.output_debug_images) {
		vector<Mat> plgs_imgs = draw_plgs(input_edge_imgs, plgs);

		// Output image representation of produced PolyLine Graphs
		write_images(emip.em_out_folder, sfm_data, plgs_imgs, "plgs_imgs_");
		for(const auto &plgs_img : plgs_imgs)
			plgs_img.release();
	}

	PolyLineGraph3DHMapImpl plg3d = PolyLineGraph3DHMapImpl();
	PLGMatchesManager plgmm(plgs,plg3d);

	vector<PolyLine2DMapSearch> plmaps;
	for(const auto &plg: plgs)
		plmaps.push_back(PolyLine2DMapSearch(plg,imgs[0].size(),4.0));

	Mat** all_fundamental_matrices = generate_all_fundamental_matrices(sfm_data);

	EdgeManager *em = new PLGEdgeManager(imgs,sfm_data,all_fundamental_matrices,plgs, DETECTION_STARTING_RADIUS, DETECTION_CORRESPONDENCES_MULTIPLICATION_FACTOR);

	vector<Mat> edge_imgs = em->get_edge_images_original_background();

	vector<Mat> &drawn_imgs = edge_imgs;

	data_bundle mfc;
	mfc.em = em;
	mfc.cm = new PLGPCM3ViewsPLGFollowing(drawn_imgs,sfm_data, img_size,all_fundamental_matrices,plgs,plmaps);
	mfc.sfmd = &sfm_data;
	mfc.original_img_size = Size(drawn_imgs[0].cols,drawn_imgs[0].rows);
	mfc.all_fundamental_matrices = all_fundamental_matrices;
	mfc.imgs = drawn_imgs;
	mfc.outfolder = emip.em_out_folder;
	mfc.plgs = plgs;

	int first_edgepoint = sfm_data.points_.size();

	// Run edge reconstruction pipeline
	edge_reconstruction_pipeline(emip,imgs, emip.em_out_folder, plgs, sfm_data,&mfc,plgmm,plmaps);

	// Output SfM data with integrated edge points resulting from 3D edges sampling during reconstruction
	output_sfm_data(emip.sfm_data_file, sfm_data, emip.em_out_folder + string("before_filtering.json"));

	// Remove outliers in current SfM data
	filter(sfm_data, first_edgepoint);

	// Output final filtered point cloud
	output_sfm_data(emip.sfm_data_file, sfm_data, emip.output_json);

	// Output points projections images
	if(emip.output_debug_images) {
		// on original images
		draw_sfmd_points(imgs, plgs, sfm_data, first_edgepoint, emip.em_out_folder, "output_on_imgs_");
		// on original images
		draw_sfmd_points_plgs(imgs, plgs, sfm_data, first_edgepoint, emip.em_out_folder, "output_on_plgs_");
	}

	return 0;
}
