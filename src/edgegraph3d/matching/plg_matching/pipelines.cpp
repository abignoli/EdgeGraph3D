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


#include "pipelines.hpp"

#include <omp.h>
#include <opencv2/core/mat.hpp>
#include <opencv2/core/mat.inl.hpp>
#include <opencv2/core/types.hpp>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "polyline_graph_2d_hmap_impl.hpp"
#include "SfMData.h"
#include "plg_matches_manager.hpp"
#include "plg_matching_from_refpoints.hpp"
#include "polyline_matching.hpp"
#include "polyline_2d_map_search.hpp"
#include "polyline_matcher.hpp"
#include "polyline_graph_2d.hpp"
#include "polyline_graph_3d.hpp"
#include "polyline_graph_3d_hmap_impl.hpp"
#include "datatypes.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "global_defines.hpp"
#include "glm.hpp"

#include "data_bundle.hpp"
#include "drawing_utilities.hpp"
#include "edge_graph_3d_input_params.hpp"
#include "filtering_close_plgps.hpp"
#include "output_utilities.hpp"

int min(const int a, const int b) {
	return a < b ? a : b;
}

void pipeline_polyline_matching_similarity_graph(const edge_matcher_input_params &emip, vector<Mat> &imgs, const char *out_folder, const vector<PolyLineGraph2DHMapImpl> &plgs, SfMData &sfmd, data_bundle *mfc, PLGMatchesManager &plgmm, const vector<PolyLine2DMapSearch> &plmaps,
		double &pmsg_time, ulong &pmsg_matches_amount, double &refpfrom_pmsg_time, ulong &refpoints_from_pmsg_amount, vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &p3ds) {
	double pmsg_time_start = omp_get_wtime();
	tuple<vector<vector<set<ulong>>>,vector<vector<vector<ulong>>>,vector<vector<set<ulong>>>> pmsg_res = polyline_matching_similarity_graph(plgs, sfmd, mfc->original_img_size, (out_folder + string(COMPATIBILITY_GRAPH)).c_str(), (out_folder + string(POLYLINE_COMMUNITIES)).c_str());
	double pmsg_time_end = omp_get_wtime();
	pmsg_time = pmsg_time_end - pmsg_time_start;

	// Run polyline matching using similarity graph and community detection
	vector<vector<set<ulong>>> pmsg_matches = get<2>(pmsg_res);

	cout << "Found " << pmsg_matches.size() << " potential matches using similarity graphs and community detection. (1/2)\n";

	pmsg_matches_amount = pmsg_matches.size();

	// Output polyline matching images
	vector<Mat> out_pmsg_imgs;
	if(emip.output_debug_images) {
		for(int i=0; i < plgs.size(); i++)
			out_pmsg_imgs.push_back(draw_polyline_graph_simplified(imgs[i],plgs[i],Scalar(255,255,255)));
		draw_polyline_matches(out_pmsg_imgs,plgs,pmsg_matches);
		write_images(out_folder, sfmd, out_pmsg_imgs, "pmsg_");
	}

	// Run PLG matching using current polyline matches
	double refpfrom_pmsg_time_start = omp_get_wtime();
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> pmsg_p3ds;
	for(int i=0; i < pmsg_matches.size(); i++) {
		const vector<set<ulong>> &potentially_compatible_polylines = pmsg_matches[i];
		cout << "Extracting 3D edges from potential polyline match " << i << "...\n";
#ifdef SWITCH_RUNPARALLEL
		vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> cur_pmsg_p3ds = find_new_3d_points_from_compatible_polylines_expandallviews_parallel(sfmd, plgs, mfc->all_fundamental_matrices, potentially_compatible_polylines,plgmm,plmaps);
#else
		vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> cur_pmsg_p3ds = find_new_3d_points_from_compatible_polylines_expandallviews(sfmd, plgs, mfc->all_fundamental_matrices, potentially_compatible_polylines,plgmm,plmaps);
#endif
		if(cur_pmsg_p3ds.size() > 0)
			pmsg_p3ds.insert(pmsg_p3ds.end(),cur_pmsg_p3ds.begin(),cur_pmsg_p3ds.end());
	}
	double refpfrom_pmsg_time_end = omp_get_wtime();
	refpfrom_pmsg_time = refpfrom_pmsg_time_end - refpfrom_pmsg_time_start;

	refpoints_from_pmsg_amount = pmsg_p3ds.size();

	p3ds = pmsg_p3ds;
}

void pipeline_polyline_matching_closeness_to_refpoints(const edge_matcher_input_params &emip, vector<Mat> &imgs, const char *out_folder, const vector<PolyLineGraph2DHMapImpl> &plgs, SfMData &sfmd, data_bundle *mfc, PLGMatchesManager &plgmm, const vector<PolyLine2DMapSearch> &plmaps,
		double &pmctr_time, ulong &pmctr_matches_amount, double &refpfrom_pmctr_time, ulong &refpoints_from_pmctr_amount, vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &p3ds) {

	// Run polyline matching using Euclidean distance from reference points
	double pmctr_time_start = omp_get_wtime();
	pair<vector<ulong>,vector<vector<set<ulong>>>> pmctr_res = polyline_matching_closeness_to_refpoints(plgs, sfmd,mfc->original_img_size);
	double pmctr_time_end = omp_get_wtime();
	pmctr_time = pmctr_time_end - pmctr_time_start;

	vector<vector<set<ulong>>> pmctr_matches = pmctr_res.second;

	pmctr_matches_amount = pmctr_matches.size();

	cout << "Found " << pmctr_matches.size() << " potential matches using Euclidean distance from reference points. (2/2)\n";

	if(emip.output_debug_images) {
		vector<Mat> out_pmctr_imgs;
		for(int i=0; i < plgs.size(); i++)
			out_pmctr_imgs.push_back(draw_polyline_graph_simplified(imgs[i],plgs[i],Scalar(255,255,255)));

		draw_polyline_matches(out_pmctr_imgs,plgs,pmctr_matches);
		write_images(out_folder, sfmd, out_pmctr_imgs, "pmctr_");
	}

	// Run PLG matching using current polyline matches
	double refpfrom_pmctr_time_start = omp_get_wtime();
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> pmctr_p3ds;
	for(int i=0; i < pmctr_matches.size(); i++) {
		const vector<set<ulong>> &potentially_compatible_polylines = pmctr_matches[i];
		cout << "Extracting 3D edges from potential polyline match " << i << "...\n";
#ifdef SWITCH_RUNPARALLEL
		vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> cur_pmctr_p3ds = find_new_3d_points_from_compatible_polylines_expandallviews_parallel(sfmd, plgs, mfc->all_fundamental_matrices, potentially_compatible_polylines,plgmm,plmaps);
#else
		vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> cur_pmctr_p3ds = find_new_3d_points_from_compatible_polylines_expandallviews(sfmd, plgs, mfc->all_fundamental_matrices, potentially_compatible_polylines,plgmm,plmaps);
#endif
		if(cur_pmctr_p3ds.size() > 0)
			pmctr_p3ds.insert(pmctr_p3ds.end(),cur_pmctr_p3ds.begin(),cur_pmctr_p3ds.end());
	}
	double refpfrom_pmctr_time_end = omp_get_wtime();
	refpfrom_pmctr_time = refpfrom_pmctr_time_end - refpfrom_pmctr_time_start;

	refpoints_from_pmctr_amount = pmctr_p3ds.size();

	for(const auto &p: pmctr_p3ds)
		p3ds.push_back(p);
}

void pipeline_refpoints(const edge_matcher_input_params &emip, vector<Mat> &imgs, const char *out_folder, const vector<PolyLineGraph2DHMapImpl> &plgs, SfMData &sfmd, data_bundle *mfc, PLGMatchesManager &plgmm, const vector<PolyLine2DMapSearch> &plmaps,
		double &refpfrom_refpoints_time, ulong &refpoints_from_refpoints_amount, vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &p3ds) {
	double refpfrom_refpoints_time_start = omp_get_wtime();
#ifdef SWITCH_RUNPARALLEL
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> p3ds_r = plg_matching_from_refpoints_parallel(sfmd, mfc->em, mfc->cm,plgmm);
#else
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> p3ds_r = plg_matching_from_refpoints(sfmd, mfc->em, mfc->cm,plgmm);
#endif

	double refpfrom_refpoints_time_end = omp_get_wtime();
	refpfrom_refpoints_time = refpfrom_refpoints_time_end - refpfrom_refpoints_time_start;

	refpoints_from_refpoints_amount = p3ds_r.size();

	for(const auto &p: p3ds_r)
		p3ds.push_back(p);
}

void print_final_stats(SfMData &sfmd, const ulong refpoints_amount, const ulong pmsg_matches_amount,
		const ulong refpoints_from_pmsg_amount, const ulong pmctr_matches_amount, const ulong refpoints_from_pmctr_amount,
		const ulong refpoints_from_refpoints_amount, const double pmsg_time, const double refpfrom_pmsg_time,
		const double pmctr_time, const double refpfrom_pmctr_time, const double refpfrom_refpoints_time) {
	cout << "******************************************" << endl;

	cout << "Initial reference points: " << refpoints_amount << endl;
	cout << "Final 3D points: " << sfmd.numPoints_ << endl;
	cout << "Polyline matches (1/2): " << pmsg_matches_amount << endl;
	cout << "3D points from polyline matches (1/2): " << refpoints_from_pmsg_amount << endl;
	cout << "Polyline matches (2/2) amount: " << pmctr_matches_amount << endl;
	cout << "3D points from polyline matches (2/2): " << refpoints_from_pmctr_amount << endl;
	cout << "3D points from reference points: " << refpoints_from_refpoints_amount << endl;
	cout << endl;
	cout << "Polyline matching (1/2): " << pmsg_time << "s" << endl;
	cout << "3D points extraction from polyline matches (1/2): " << refpfrom_pmsg_time << "s" << endl;
	cout << "Polyline matching (1/2): " << pmctr_time << "s" << endl;
	cout << "3D points extraction from polyline matches (2/2): " << refpfrom_pmctr_time << "s" << endl;
	cout << "3D points extraction from reference points: " << refpfrom_refpoints_time << "s" << endl;

	cout << "******************************************" << endl;
}

void edge_reconstruction_pipeline(const edge_matcher_input_params &emip, vector<Mat> &imgs, const char *out_folder, const vector<PolyLineGraph2DHMapImpl> &plgs, SfMData &sfmd, data_bundle *mfc, PLGMatchesManager &plgmm, const vector<PolyLine2DMapSearch> &plmaps)
{
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> p3ds;

	double pmsg_time, refpfrom_pmsg_time;
	ulong pmsg_matches_amount, refpoints_from_pmsg_amount;

	double pmctr_time, refpfrom_pmctr_time;
	ulong pmctr_matches_amount, refpoints_from_pmctr_amount;

	double refpfrom_refpoints_time;
	ulong refpoints_from_refpoints_amount;

	double DO_MEASURE_OMP_TIME_START_start = omp_get_wtime();

	ulong refpoints_amount = sfmd.numPoints_;

	// Run pipeline for edge reconstruction using polyline matching using similarity graph and community detection
	pipeline_polyline_matching_similarity_graph(emip, imgs, out_folder, plgs, sfmd, mfc, plgmm, plmaps,
			pmsg_time, pmsg_matches_amount, refpfrom_pmsg_time, refpoints_from_pmsg_amount, p3ds);

	// Run pipeline for edge reconstruction using polyline matching using Euclidean distance from reference points
	pipeline_polyline_matching_closeness_to_refpoints(emip, imgs, out_folder, plgs, sfmd, mfc, plgmm, plmaps,
			pmctr_time, pmctr_matches_amount, refpfrom_pmctr_time, refpoints_from_pmctr_amount, p3ds);

	// Run pipeline for edge reconstruction using reference points
	pipeline_refpoints(emip, imgs, out_folder, plgs, sfmd, mfc, plgmm, plmaps,
			refpfrom_refpoints_time, refpoints_from_refpoints_amount, p3ds);

	double DO_MEASURE_OMP_TIME_END_end = omp_get_wtime();

	// Serialize and save PLG matches manager, including actual 3D edges structure
	serialize_plg(plgmm.get_plg3d(), out_folder + string(PLG3D_OUTNAME));

	// Limit density of added edge-points to a proper extent
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> filtered_p3ds = filter_3d_points_close_2d_array(plgs,imgs[0].size(),p3ds);

	// Add computed edge-points to original SfM data
	add_3dpoints_to_sfmd(sfmd,filtered_p3ds);

	// Print final statistics
	print_final_stats(sfmd, refpoints_amount, pmsg_matches_amount,
			refpoints_from_pmsg_amount, pmctr_matches_amount, refpoints_from_pmctr_amount,
			refpoints_from_refpoints_amount, pmsg_time, refpfrom_pmsg_time,
			pmctr_time, refpfrom_pmctr_time, refpfrom_refpoints_time);
}


