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


#include "polyline_matcher.hpp"

#include <stddef.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>

#include "SfMData.h"
#include "polyline_2d_map_search.hpp"
#include "community_detection_interface.hpp"
#include "graph_adjacency_set_undirected_no_type_weighted.hpp"
#include "global_defines.hpp"
#include "glm.hpp"



struct KeyFuncs_cam_polyline
{
    size_t operator()(const pair<int,ulong>& k)const
    {
        return std::hash<int>()(k.first) ^ std::hash<ulong>()(k.second);
    }

    bool operator()(const pair<int,ulong>& k1, const pair<int,ulong>& k2)const
    {
            return k1.first == k2.first && k1.second == k2.second;
    }
};

ulong polyline_matching4_get_polyline_graph_node(GraphAdjacencySetUndirectedNoType &pmg, unordered_map<pair<int,ulong>,ulong,KeyFuncs_cam_polyline,KeyFuncs_cam_polyline> &polyline_matches_map, vector<pair<int,ulong>> &polyline_matches_vector, const pair<int,ulong> &pl) {
	unordered_map<pair<int,ulong>,ulong,KeyFuncs_cam_polyline,KeyFuncs_cam_polyline>::iterator it = polyline_matches_map.find(pl);
	ulong node_id;
	if(it != polyline_matches_map.end())
		node_id = it->second;
	else {
		node_id = pmg.add_node();
		polyline_matches_map[pl] = node_id;
		polyline_matches_vector.push_back(pl);
	}
	return node_id;
}

pair<vector<ulong>,vector<vector<set<ulong>>>> polyline_matching_closeness_to_refpoints(const vector<PolyLineGraph2DHMapImpl> &plgs, const SfMData &sfmd, const Size &img_sz)
{
	cout << "Polyline matching (2/2)...\n";

	vector<PolyLine2DMapSearch> plmaps;
	for(const auto &plg: plgs)
		plmaps.push_back(PolyLine2DMapSearch(plg,img_sz,FIND_WITHIN_DIST));

	DO_MEASURE_TIME_START

	bool tmp;
	vector<ulong> refpoints;

	ulong pl_a,pl_b;

	unordered_map<pair<int,ulong>,ulong,KeyFuncs_cam_polyline,KeyFuncs_cam_polyline> polyline_matches_map;
	vector<pair<int,ulong>> polyline_matches_vector;

	GraphAdjacencySetUndirectedNoType pmg;

	for(ulong refpoint_id = 0; refpoint_id < sfmd.numPoints_ ; refpoint_id++)
	{
		// cout << "Processing refpoint " << refpoint_id << "\n";

		vector<vector<tuple<ulong,ulong,glm::vec2,float>>> curpoint_res;

		for(const auto cam_id : sfmd.camViewingPointN_[refpoint_id]) {
			const PolyLineGraph2DHMapImpl &plg = plgs[cam_id];
			curpoint_res.push_back(plmaps[cam_id].find_polylines_within_search_dist_with_reprojections(get_2d_coordinates_of_point_on_image(sfmd,cam_id,refpoint_id,tmp)));
		}

		int maxpl = 0;
		for(const auto &cpr : curpoint_res)
			maxpl = maxpl < cpr.size() ? cpr.size() : maxpl;

		if(maxpl == 1) {
			set<pair<int,ulong>> cur_cams_pls;

			float min_dist = std::numeric_limits<float>::max();
			float max_dist = std::numeric_limits<float>::min();

			for(int i=0; i < curpoint_res.size(); i++) {
				if(curpoint_res[i].size() == 0)
					continue;

				const tuple<ulong,ulong,glm::vec2,float> &pl = curpoint_res[i][0];
				min_dist = min_dist <= get<3>(pl) ? min_dist : get<3>(pl);
				max_dist = max_dist >= get<3>(pl) ? max_dist : get<3>(pl);
				cur_cams_pls.insert(make_pair(sfmd.camViewingPointN_[refpoint_id][i],get<0>(pl)));
			}

			if(cur_cams_pls.size() < sfmd.camViewingPointN_[refpoint_id].size() * 0.7)
				continue;

			if(min_dist < (max_dist / DETECTION_CORRESPONDENCES_MULTIPLICATION_FACTOR))
				continue;

			if(max_dist > (min_dist * DETECTION_CORRESPONDENCES_MULTIPLICATION_FACTOR))
				continue;

			if(cur_cams_pls.size() < 2)
				continue;

			vector<ulong> pl_ids;
			for(set<pair<int,ulong>>::iterator i=cur_cams_pls.begin(); i != cur_cams_pls.end();i++)
				pl_ids.push_back(polyline_matching4_get_polyline_graph_node(pmg,polyline_matches_map,polyline_matches_vector,*i));

			for(int i=0; i < cur_cams_pls.size(); i++)
				for(int j=i+1; j < cur_cams_pls.size(); j++)
					pmg.add_edge(pl_ids[i],pl_ids[j]);



			refpoints.push_back(refpoint_id);
		}
	}

	vector<vector<ulong>> components = pmg.get_components();
	vector<vector<set<ulong>>> polyline_matches;
	for(const auto &component : components) {
		vector<set<ulong>> cur_component_pls(plgs.size());
		for(const auto pmg_node : component) {
			const pair<int,ulong> &cur_pl = polyline_matches_vector[pmg_node];
			cur_component_pls[cur_pl.first].insert(cur_pl.second);
		}
		polyline_matches.push_back(cur_component_pls);
	}

	// cout << refpoints.size() << " successful refpoints\n";

	DO_MEASURE_TIME_END

	return make_pair(refpoints,polyline_matches);
}

// Assumes close_refpoints_a and close_refpoints_b are both sorted
float compute_compatibility(const vector<float> &refpoints_weights, const vector<ulong> &close_refpoints_a_visible_on_img_b, const vector<ulong> &close_refpoints_b_visible_on_img_a) {
	std::vector<ulong>::iterator it,it2;
	vector<ulong> v(close_refpoints_a_visible_on_img_b.size()+close_refpoints_b_visible_on_img_a.size());

	it=std::set_intersection (close_refpoints_a_visible_on_img_b.begin(), close_refpoints_a_visible_on_img_b.end(),close_refpoints_b_visible_on_img_a.begin(), close_refpoints_b_visible_on_img_a.end(),v.begin());
	float intersection_weightsum=0.0;
	for(it2=v.begin(); it2 != it; it2++)
		intersection_weightsum += refpoints_weights[*it2];

	if(intersection_weightsum == 0.0)
		return 0.0;

	it=std::set_union (close_refpoints_a_visible_on_img_b.begin(), close_refpoints_a_visible_on_img_b.end(),close_refpoints_b_visible_on_img_a.begin(), close_refpoints_b_visible_on_img_a.end(),v.begin());
	float union_weightsum=0.0;
	for(it2=v.begin(); it2 != it; it2++)
		union_weightsum += refpoints_weights[*it2];

	return intersection_weightsum / union_weightsum;
}

float compute_refpoint_weight(const vector<set<ulong>> &close_polylines) {
	int non_empty=0;
	int sum_pls=0;
	for(const auto &spls : close_polylines)
		if(spls.size() > 0) {
			non_empty++;
			sum_pls += spls.size();
		}
	return non_empty == 0 ? 0.0 : non_empty / ((float) sum_pls);
}

vector<vector<set<ulong>>> compute_polyline_matches_from_nodes_component_ids(const vector<pair<int,ulong>> &polyline_matches_vector,int amount_of_plgs,const vector<long> &nodes_component_ids) {
	ulong num_components = vec_max(nodes_component_ids)+1;

	vector<vector<set<ulong>>> res;
	for(ulong i=0; i < num_components; i++)
		res.push_back(vector<set<ulong>>(amount_of_plgs));
	for(ulong i=0; i < polyline_matches_vector.size(); i++)
		if(nodes_component_ids[i] >= 0) {
			const pair<int,ulong> &m = polyline_matches_vector[i];
			res[nodes_component_ids[i]][m.first].insert(m.second);
		}
	return res;
}

/* Compute close components:
 * Output:
 *	- set of polyline_id for each close polyline to refpoint on img, for each img, for each refpoint
 *	- vector of close refpoints (sorted in ascending order), for each polyline, for each image
 *	- for each component in the potential compatibility graph (i.e. elements both close to at least a refpoint): for each image the potential matches (polyline ids)
 */
tuple<vector<vector<set<ulong>>>,vector<vector<vector<ulong>>>,vector<vector<set<ulong>>>> polyline_matching_similarity_graph(const vector<PolyLineGraph2DHMapImpl> &plgs, const SfMData &sfmd, const Size &img_sz, const char *out_compatibility_graph_file, const char *out_polyline_communities)
{
	vector<PolyLine2DMapSearch> plmaps;
	for(const auto &plg: plgs)
		plmaps.push_back(PolyLine2DMapSearch(plg,img_sz,FIND_WITHIN_DIST));

	cout << "Matching 2D polylines...\n";

	DO_MEASURE_TIME_START

	vector<vector<set<ulong>>> close_polylines;
	vector<vector<vector<ulong>>> close_refpoints;
	for(int i=0; i < plgs.size(); i++)
		close_refpoints.push_back(vector<vector<ulong>>(plgs[i].polylines.size()));

	bool tmp;

	unordered_map<pair<int,ulong>,ulong,KeyFuncs_cam_polyline,KeyFuncs_cam_polyline> polyline_matches_map;
	vector<pair<int,ulong>> polyline_matches_vector;

	GraphAdjacencySetUndirectedNoType pmg;

	for(ulong refpoint_id = 0; refpoint_id < sfmd.numPoints_ ; refpoint_id++)
	{
		vector<set<ulong>> cur_close_polylines(plgs.size());

		vector<vector<tuple<ulong,ulong,glm::vec2,float>>> curpoint_res;

		for(const auto cam_id : sfmd.camViewingPointN_[refpoint_id]) {
			const PolyLineGraph2DHMapImpl &plg = plgs[cam_id];
			curpoint_res.push_back(plmaps[cam_id].find_polylines_within_search_dist_with_reprojections(get_2d_coordinates_of_point_on_image(sfmd,cam_id,refpoint_id,tmp)));
		}

		set<pair<int,ulong>> cur_cams_pls;
		for(int i=0; i < curpoint_res.size(); i++) {
			if(curpoint_res[i].size() == 0)
				continue;

			for(const auto &pl : curpoint_res[i])
				cur_cams_pls.insert(make_pair(sfmd.camViewingPointN_[refpoint_id][i],get<0>(pl)));
		}

		for(const auto &cc : cur_cams_pls) {
			cur_close_polylines[cc.first].insert(cc.second);
			close_refpoints[cc.first][cc.second].push_back(refpoint_id);
		}

		vector<ulong> pl_ids;
		for(set<pair<int,ulong>>::iterator i=cur_cams_pls.begin(); i != cur_cams_pls.end();i++)
			pl_ids.push_back(polyline_matching4_get_polyline_graph_node(pmg,polyline_matches_map,polyline_matches_vector,*i));

		for(int i=0; i < cur_cams_pls.size(); i++)
			for(int j=i+1; j < cur_cams_pls.size(); j++)
				pmg.add_edge(pl_ids[i],pl_ids[j]);

		close_polylines.push_back(cur_close_polylines);
	}


	vector<set<ulong>> point_sets_visible_from_camN;
	for(int i=0; i < plgs.size(); i++) {
		set<ulong> cs;

		for(const auto p: sfmd.pointsVisibleFromCamN_[i])
			cs.insert(p);

		point_sets_visible_from_camN.push_back(cs);
	}

	/**
	 * For img_1 in imgs:
	 * 		for pl in img1.pls:
	 * 			for img2 in imgs:
	 * 				vector of refpoints visible on img1 and img2, close to pl on img1 (sorded in ascending order)
	 */
	vector<vector<vector<vector<ulong>>>> close_refpoints_divided_by_visibility;
	for(int img1=0; img1 < plgs.size(); img1++) {
		close_refpoints_divided_by_visibility.push_back(vector<vector<vector<ulong>>>(plgs[img1].polylines.size()));
		for(int pl_id=0; pl_id < plgs[img1].polylines.size(); pl_id++) {
			close_refpoints_divided_by_visibility[img1][pl_id]=vector<vector<ulong>>(plgs.size());
			for(const auto refpoint_id : close_refpoints[img1][pl_id])
				for(int img2=0; img2 < plgs.size(); img2++)
					if(is_in(point_sets_visible_from_camN[img2],refpoint_id))
						close_refpoints_divided_by_visibility[img1][pl_id][img2].push_back(refpoint_id);
		}
	}

	// Compute refpoints' weights
	vector<float> refpoints_weights(sfmd.numPoints_);
	for(ulong refpoint_id = 0; refpoint_id < sfmd.numPoints_; refpoint_id++)
		refpoints_weights[refpoint_id] = compute_refpoint_weight(close_polylines[refpoint_id]);

	GraphAdjacencySetUndirectedNoTypeWeighted pmgw(pmg.get_nodes_num());
	float w;
	for(ulong node1=0; node1 < pmg.get_nodes_num(); node1++)
		for(const auto node2: pmg.adjacency_lists[node1])
			if(node1 < node2) {
				// Compute compatibility for node1, node2
				const pair<int,ulong> &match1 = polyline_matches_vector[node1];
				const pair<int,ulong> &match2 = polyline_matches_vector[node2];

				w = compute_compatibility(refpoints_weights, close_refpoints_divided_by_visibility[match1.first][match1.second][match2.first], close_refpoints_divided_by_visibility[match2.first][match2.second][match1.first]);

				if(w > 0.0)
					pmgw.add_edge(node1,node2,w);
			}

	// Compute communities
	vector<long> nodes_component_ids = compute_communities(pmgw, out_compatibility_graph_file, out_polyline_communities);

	// Extract matches
	vector<vector<set<ulong>>> polyline_matches = compute_polyline_matches_from_nodes_component_ids(polyline_matches_vector,plgs.size(), nodes_component_ids);

	return make_tuple(close_polylines,close_refpoints,polyline_matches);
}
