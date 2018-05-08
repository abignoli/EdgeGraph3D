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

#include "convert_edge_images_pixel_to_segment.hpp"

#include <opencv2/core/mat.hpp>
#include <opencv2/core/mat.inl.hpp>
#include <opencv2/core/matx.hpp>
#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include <glm.hpp>
#include "polyline_graph_2d_hmap_impl.hpp"
#include "graph.hpp"
#include "graph_adjacency_list_undirected.hpp"
#include "graph_adjacency_set_undirected_no_type.hpp"
#include "polyline_graph_2d.hpp"
#include "drawing_utilities.hpp"
#include "edge_graph_3d_utilities.hpp"


using namespace std;
using namespace cv;

/**
 * Takes an edge image, assigns a Node ID to each edge pixel in the image
 */
ulong** convertEdgeImagesPixelToNodes(const Mat &img, const Vec<unsigned char, 3> &edge_color, ulong &amount_of_nodes, vector<glm::vec2> &node_coords) {
	ulong** node_id_map = create_2D_array<ulong>(img.rows, img.cols);

	node_coords.clear();

	ulong cur_edgepixel_id = 0;
	for(ulong i=0; i < img.rows; i++)
		for(ulong j=0; j < img.cols; j++)
			if(img.at<Vec3b>(i,j) == edge_color) {
				node_id_map[i][j] = cur_edgepixel_id;
				node_coords.push_back(glm::vec2(j+0.5,i+0.5));
				cur_edgepixel_id++;
			}

	amount_of_nodes = cur_edgepixel_id;

	return node_id_map;
}

GraphAdjacencyListUndirected<glm::vec2> convertEdgeImagePixelToGraph(const Mat &img, const Vec<unsigned char, 3> &edge_color) {
	ulong amount_of_nodes;
	vector<glm::vec2> node_coords;
	ulong** node_id_map = convertEdgeImagesPixelToNodes(img, edge_color,amount_of_nodes,node_coords);
	int cx,cy;
	int cur_P_node_id, cur_C_node_id;

	GraphAdjacencyListUndirected<glm::vec2> edges_graph(amount_of_nodes);
	for(ulong i = 0; i < amount_of_nodes; i++)
		edges_graph.set_node_data(i,node_coords[i]);

	for(int i=0; i < img.rows - 1; i++)
		for(int j=0; j < img.cols - 1; j++)
			if(img.at<Vec3b>(i,j) == edge_color) {
				/**
				 * Current pixel: P
				 * Pixels to check: C
				 *
				 * P  C1
				 * C2 C3
				 *
				 * if P is edge pixel and C is edge pixel, add edge P <-> C
				 */

				cur_P_node_id = node_id_map[i][j];

				// C1

				cx = j+1;
				cy = i;
				if(img.at<Vec3b>(cy,cx) == edge_color) {
					cur_C_node_id = node_id_map[cy][cx];
					edges_graph.add_edge(cur_P_node_id,cur_C_node_id);
				}

				// C2

				cx = j;
				cy = i+1;
				if(img.at<Vec3b>(cy,cx) == edge_color) {
					cur_C_node_id = node_id_map[cy][cx];
					edges_graph.add_edge(cur_P_node_id,cur_C_node_id);
				}

				// C3

				cx = j+1;
				cy = i+1;
				if(img.at<Vec3b>(cy,cx) == edge_color) {
					cur_C_node_id = node_id_map[cy][cx];
					edges_graph.add_edge(cur_P_node_id,cur_C_node_id);
				}

				if(j > 1) {
					/**
					 * X P
					 * C x
					 */

					cx = j-1;
					cy = i+1;
					if(img.at<Vec3b>(cy,cx) == edge_color) {
						cur_C_node_id = node_id_map[cy][cx];
						edges_graph.add_edge(cur_P_node_id,cur_C_node_id);
					}
				}

			}

	delete_2D_array(node_id_map);

	cout << "Graph nodes: " << edges_graph.get_nodes_num() << endl;

	return edges_graph;
}

vector<glm::vec4> convert_graph_to_edge_list(const GraphAdjacencyListUndirected<glm::vec2> &g) {
	vector<vector<ulong>> adjacency_list = g.get_adjacency_lists();
	vector<glm::vec2> nodes_data = g.get_nodes_data();
	vector<glm::vec4> res;
	ulong nodes_amount = g.get_nodes_num();

	for(ulong start_node=0;start_node<nodes_amount;start_node++) {
		const glm::vec2 start_coords = nodes_data[start_node];
		for(const auto end_node : adjacency_list[start_node])
			if(start_node <= end_node) {
				const glm::vec2 end_coords = nodes_data[end_node];
				res.push_back(glm::vec4(start_coords[0],start_coords[1],end_coords[0],end_coords[1]));
			}
	}

	return res;
}

vector<glm::vec4> convertEdgeImagePixelToSegments(const Mat &img, const Vec<unsigned char, 3> &edge_color) {
	GraphAdjacencyListUndirected<glm::vec2> g = convertEdgeImagePixelToGraph(img,edge_color);
	return convert_graph_to_edge_list(g);
}

vector<vector<glm::vec4>> convertEdgeImagePixelToSegments(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color) {
	vector<vector<glm::vec4>> res;
	for(const auto &img: imgs)
		res.push_back(convertEdgeImagePixelToSegments(img,edge_color));
	return res;
}

Mat convertEdgeImagePixelToSegmentsImage(const Mat &img, const Vec<unsigned char, 3> &edge_color) {
	vector<glm::vec4> segments = convertEdgeImagePixelToSegments(img,edge_color);
	cout << "Number of segments: " << segments.size() << endl;
	Mat outimg = get_black_image(img);
	draw_segments_on_image(outimg,segments);
	return outimg;
}

vector<Mat> convertEdgeImagesPixelToSegmentsImages(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color) {
	vector<Mat> res;
		for(const auto &img: imgs)
			res.push_back(convertEdgeImagePixelToSegmentsImage(img,edge_color));
	return res;
}

vector<GraphAdjacencyListUndirected<glm::vec2>> convertEdgeImagesPixelToGraphs(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color) {
	vector<GraphAdjacencyListUndirected<glm::vec2>> res;
	for(const auto &img: imgs)
		res.push_back(convertEdgeImagePixelToGraph(img,edge_color));
	return res;
}

inline bool is_edge(const Mat &img,const int i, const int j, const Vec<unsigned char, 3> &edge_color) {
	return img.at<Vec3b>(i,j) == edge_color;
}

/**
 * Takes an edge image, assigns a Node ID to each edge pixel in the image
 */
ulong** convertEdgeImagesPixelToNodesNoSquaresNoTriangles(const Mat &img, const Vec<unsigned char, 3> &edge_color, ulong &amount_of_nodes, vector<glm::vec2> &node_coords) {
	ulong** node_id_map = create_2D_array<ulong>(img.rows, img.cols);
	bool** set_node_id_map = create_2D_array<bool>(img.rows, img.cols);

	node_coords.clear();

	for(int i=0; i < img.rows; i++)
		for(int j=0; j < img.cols; j++)
			set_node_id_map[i][j] = false;

	ulong cur_edgepixel_id = 0;
	for(int i=0; i < img.rows; i++)
		for(int j=0; j < img.cols; j++)
			if(is_edge(img,i,j,edge_color) && !set_node_id_map[i][j]) {
				glm::vec2 new_coords;
				if((i < img.rows -1) && (j < img.cols -1) && (is_edge(img,i+1,j,edge_color) && is_edge(img,i,j+1,edge_color) && is_edge(img,i+1,j+1,edge_color))) // Check for squares
				{
					/*
					 *   P *
					 *   * *
					 */
					new_coords = glm::vec2(j+1,i+1);
					set_node_id_map[i][j+1]=true;
					set_node_id_map[i+1][j]=true;
					set_node_id_map[i+1][j+1]=true;
					node_id_map[i][j+1] = cur_edgepixel_id;
					node_id_map[i+1][j] = cur_edgepixel_id;
					node_id_map[i+1][j+1] = cur_edgepixel_id;
				} else if ((i < img.rows -1) && (j < img.cols -1) && (is_edge(img,i,j+1,edge_color) && is_edge(img,i+1,j+1,edge_color))) // Check for triangle 1
				{
					/*
					 *   P *
					 *     *
					 */
					new_coords = glm::vec2(j+0.5+2.0/3,i+0.5+1.0/3);
					set_node_id_map[i][j+1]=true;
					set_node_id_map[i+1][j+1]=true;
					node_id_map[i][j+1] = cur_edgepixel_id;
					node_id_map[i+1][j+1] = cur_edgepixel_id;

				} else if ((i < img.rows -1) && (j < img.cols -1) && (is_edge(img,i+1,j,edge_color) && is_edge(img,i+1,j+1,edge_color))) // Check for triangle 2
				{
					/*
					 *   P
					 *   * *
					 */
					new_coords = glm::vec2(j+0.5+1.0/3,i+0.5+2.0/3);
					set_node_id_map[i+1][j]=true;
					set_node_id_map[i+1][j+1]=true;
					node_id_map[i+1][j] = cur_edgepixel_id;
					node_id_map[i+1][j+1] = cur_edgepixel_id;

				} else if ((i < img.rows -1) && (j > 0) && (is_edge(img,i+1,j-1,edge_color) && is_edge(img,i+1,j,edge_color))) // Check for triangle 2
				{
					/*
					 *     P
					 *   * *
					 */
					new_coords = glm::vec2(j+0.5-1.0/3,i+0.5+2.0/3);
					set_node_id_map[i+1][j-1]=true;
					set_node_id_map[i+1][j]=true;
					node_id_map[i+1][j-1] = cur_edgepixel_id;
					node_id_map[i+1][j] = cur_edgepixel_id;

				} else {
					new_coords = glm::vec2(j+0.5,i+0.5);
				}
				node_id_map[i][j] = cur_edgepixel_id;
				node_coords.push_back(new_coords);
				set_node_id_map[i][j]=true;
				cur_edgepixel_id++;
			}

	amount_of_nodes = cur_edgepixel_id;
	delete_2D_array<bool>(set_node_id_map);

	return node_id_map;
}

/**
 * Takes an edge image, assigns a Node ID to each edge pixel in the image
 */
ulong** convertEdgeImagesPixelToNodesNoSquaresNoTriangles_remove_useless_hubs(const Mat &c_img, const Vec<unsigned char, 3> &edge_color, ulong &amount_of_nodes, vector<glm::vec2> &node_coords) {
	Mat img = Mat(c_img);
	ulong** node_id_map = create_2D_array<ulong>(img.rows, img.cols);
	bool** set_node_id_map = create_2D_array<bool>(img.rows, img.cols);
	const Vec<unsigned char, 3> non_edge_color(edge_color[0]-1,edge_color[1]-1,edge_color[2]-1);

	node_coords.clear();

	for(int i=0; i < img.rows; i++)
		for(int j=0; j < img.cols; j++)
			set_node_id_map[i][j] = false;

	ulong cur_edgepixel_id = 0;
	for(int i=0; i < img.rows; i++)
		for(int j=0; j < img.cols; j++)
			if(is_edge(img,i,j,edge_color) && !set_node_id_map[i][j]) {
				glm::vec2 new_coords;

				/**
				 * TL TC TR
				 * LC PP RC
				 * BL BC BR
				 *
				 * disable pp if not useful
				 */
				if(
						(i>1 && j>1 && is_edge(img,i-1,j,edge_color) && is_edge(img,i,j-1,edge_color) && !is_edge(img,i+1,j+1,edge_color)) ||
						(i>1 && j<img.cols-1 && is_edge(img,i-1,j,edge_color) && is_edge(img,i,j+1,edge_color) && !is_edge(img,i+1,j-1,edge_color)) ||
						(i<img.rows-1 && j<img.cols-1 && is_edge(img,i+1,j,edge_color) && is_edge(img,i,j+1,edge_color) && !is_edge(img,i-1,j-1,edge_color)) ||
						(i<img.rows-1 && j>1 && is_edge(img,i+1,j,edge_color) && is_edge(img,i,j-1,edge_color) && !is_edge(img,i-1,j+1,edge_color))
				) {
					// Clear pixel
					img.at<Vec3b>(i,j) = non_edge_color;
				} else
				{
					new_coords = glm::vec2(j+0.5,i+0.5);
					node_id_map[i][j] = cur_edgepixel_id;
					node_coords.push_back(new_coords);
					set_node_id_map[i][j]=true;
					cur_edgepixel_id++;
				}
			}

	amount_of_nodes = cur_edgepixel_id;
	delete_2D_array<bool>(set_node_id_map);

	img.release();

	return node_id_map;
}

#define LOOP_CHECK_DIST 8

pair<GraphAdjacencySetUndirectedNoType,vector<glm::vec2>> convertEdgeImagePixelToGraph_NoCycles(const Mat &img, const Vec<unsigned char, 3> &edge_color) {
	ulong amount_of_nodes;
	vector<glm::vec2> node_coords;

	// DIsable nocycle
	//ulong** node_id_map = convertEdgeImagesPixelToNodesNoSquaresNoTriangles(img, edge_color,amount_of_nodes,node_coords);
	//ulong** node_id_map = convertEdgeImagesPixelToNodes(img, edge_color,amount_of_nodes,node_coords);
	ulong** node_id_map = convertEdgeImagesPixelToNodesNoSquaresNoTriangles_remove_useless_hubs(img, edge_color,amount_of_nodes,node_coords);

	int cx,cy;
	int cur_P_node_id, cur_C_node_id;

	GraphAdjacencySetUndirectedNoType edges_graph(amount_of_nodes);

	for(int i=0; i < img.rows - 1; i++)
		for(int j=0; j < img.cols - 1; j++)
			if(img.at<Vec3b>(i,j) == edge_color) {
				/**
				 * Current pixel: P
				 * Pixels to check: C
				 *
				 * P  C1
				 * C2 C3
				 *
				 * if P is edge pixel and C is edge pixel, add edge P <-> C
				 */

				cur_P_node_id = node_id_map[i][j];

				// C1

				cx = j+1;
				cy = i;
				if(img.at<Vec3b>(cy,cx) == edge_color) {
					cur_C_node_id = node_id_map[cy][cx];
					if(cur_P_node_id != cur_C_node_id && !edges_graph.is_connected(cur_P_node_id, cur_C_node_id, LOOP_CHECK_DIST))
						edges_graph.add_edge(cur_P_node_id,cur_C_node_id);
				}

				// C2

				cx = j;
				cy = i+1;
				if(img.at<Vec3b>(cy,cx) == edge_color) {
					cur_C_node_id = node_id_map[cy][cx];
					if(cur_P_node_id != cur_C_node_id && !edges_graph.is_connected(cur_P_node_id, cur_C_node_id, LOOP_CHECK_DIST))
						edges_graph.add_edge(cur_P_node_id,cur_C_node_id);
				}

				// C3

				cx = j+1;
				cy = i+1;
				if(img.at<Vec3b>(cy,cx) == edge_color) {
					cur_C_node_id = node_id_map[cy][cx];
					if(cur_P_node_id != cur_C_node_id && !edges_graph.is_connected(cur_P_node_id, cur_C_node_id, LOOP_CHECK_DIST))
						edges_graph.add_edge(cur_P_node_id,cur_C_node_id);
				}

				if(j > 1) {
					/**
					 * X P
					 * C x
					 */

					cx = j-1;
					cy = i+1;
					if(img.at<Vec3b>(cy,cx) == edge_color) {
						cur_C_node_id = node_id_map[cy][cx];
						if(cur_P_node_id != cur_C_node_id && !edges_graph.is_connected(cur_P_node_id, cur_C_node_id, LOOP_CHECK_DIST))
							edges_graph.add_edge(cur_P_node_id,cur_C_node_id);
					}
				}

			}

	delete_2D_array(node_id_map);

	return make_pair(edges_graph,node_coords);
}

inline bool is_isolated(const int num_adj) {
	return num_adj == 0;
}

inline bool is_ppline(const int num_adj) {
	return num_adj == 2;
}

inline bool is_extreme(const int num_adj) {
	return num_adj == 1;
}

inline bool is_hub(const int num_adj) {
	return num_adj > 2;
}

inline bool is_polylineend(const int num_adj) {
	return !is_ppline(num_adj);
}

inline bool is_isolated(const vector<std::set<ulong>> &adjacencies, const int node_id) {
	return is_isolated(adjacencies[node_id].size());
}

inline bool is_ppline(const vector<std::set<ulong>> &adjacencies, const int node_id) {
	return is_ppline(adjacencies[node_id].size());
}

inline bool is_extreme(const vector<std::set<ulong>> &adjacencies, const int node_id) {
	return is_extreme(adjacencies[node_id].size());
}

inline bool is_hub(const vector<std::set<ulong>> &adjacencies, const int node_id) {
	return is_hub(adjacencies[node_id].size());
}

inline bool is_polylineend(const vector<std::set<ulong>> &adjacencies, const int node_id) {
	return is_polylineend(adjacencies[node_id].size());
}

inline pair<ulong,ulong> get_ppline_neighbors(const std::set<ulong> &cur_adj) {
	std::set<ulong>::iterator it = cur_adj.begin();
	ulong prev = *it;
	it++;
	ulong next = *it;
	return make_pair(prev,next);
}

inline ulong get_ppline_neighbor_no_come_back(const std::set<ulong> &cur_adj, const int no_come_back_id) {
	std::set<ulong>::iterator it = cur_adj.begin();
	ulong prev = *it;
	it++;
	ulong next = *it;
	return prev != no_come_back_id ? prev : next;
}

/**
 * Find extreme without coming back to no_come_back_id
 */
void find_polylineend_no_come_back(const int start_node_id,const vector<std::set<ulong>> &adjacencies, const int no_come_back_id, vector<ulong> &res) {
	ulong cur_node_id, prev_node_id,next_node_id;

	prev_node_id = no_come_back_id;
	cur_node_id = start_node_id;
	res.push_back(cur_node_id);

	// if cur_node_id == no_come_back_id -> loop detected
	while(cur_node_id != no_come_back_id && is_ppline(adjacencies[cur_node_id].size())) {
		next_node_id = get_ppline_neighbor_no_come_back(adjacencies[cur_node_id],prev_node_id);
		prev_node_id = cur_node_id;
		cur_node_id = next_node_id;
		res.push_back(cur_node_id);
	}
}

vector<vector<ulong>> find_polylines_from_hub(const int start_node_id,const vector<std::set<ulong>> &adjacencies) {
	vector<vector<ulong>> res;
	for(std::set<ulong>::iterator it = adjacencies[start_node_id].begin(); it != adjacencies[start_node_id].end(); it++) {
		vector<ulong> cur_res;
		cur_res.push_back(start_node_id);
		find_polylineend_no_come_back(*it,adjacencies, start_node_id, cur_res);
		res.push_back(cur_res);
	}
	return res;
}

vector<ulong> find_polyline_from_ppline(const int start_node_id,const vector<std::set<ulong>> &adjacencies) {
	vector<ulong> res;
	int num_adj;
	int cur_node_id = start_node_id;
	const std::set<ulong> &cur_adj = adjacencies[cur_node_id];

	pair<ulong,ulong> ns = get_ppline_neighbors(cur_adj);
	int prev = ns.first;
	int next = ns.second;

	// find an extreme from prev
	find_polylineend_no_come_back(prev,adjacencies,cur_node_id,res);

	// reverse and add current point
	std::reverse(res.begin(), res.end());
	res.push_back(cur_node_id);

	// Loop of pplines
	if(res[0] == res[res.size()-1])
		return res;

	// find other extreme from next
	find_polylineend_no_come_back(next,adjacencies,cur_node_id,res);

	return res;
}

vector<ulong> find_polyline_from_extreme(const int start_node_id,const vector<std::set<ulong>> &adjacencies) {
	vector<ulong> res;
	int num_adj;
	ulong cur_node_id, prev_node_id;

	cur_node_id=start_node_id;
	res.push_back(cur_node_id);
	const std::set<ulong> &cur_adj = adjacencies[cur_node_id];

	prev_node_id = start_node_id;
	cur_node_id = *cur_adj.begin();
	find_polylineend_no_come_back(cur_node_id,adjacencies,prev_node_id,res);

	return res;
}

vector<vector<ulong>> find_polylines(const int start_node_id,const vector<std::set<ulong>> &adjacencies) {
	vector<vector<ulong>> res;
	const std::set<ulong> &cur_adj = adjacencies[start_node_id];
	const int num_adj = cur_adj.size();

	if(is_ppline(num_adj)) {
		// Node i is on a polyline
		res.push_back(find_polyline_from_ppline(start_node_id,adjacencies));
	} else if(is_extreme(num_adj)) {
		// Node i is an extreme
		res.push_back(find_polyline_from_extreme(start_node_id,adjacencies));
	} else if(is_hub(num_adj)) {
		// Node i is a hub
		res = find_polylines_from_hub(start_node_id, adjacencies);
	}

	return res;
}

vector<glm::vec2> get_polyline_coords(const vector<ulong> &polyline, const vector<glm::vec2> &all_coords) {
	vector<glm::vec2> polyline_coords;
	for(const auto polyline_node : polyline)
		polyline_coords.push_back(all_coords[polyline_node]);
	return polyline_coords;
}

PolyLineGraph2DHMapImpl convert_EdgeGraph_to_PolyLineGraph(const pair<GraphAdjacencySetUndirectedNoType,vector<glm::vec2>> &p) {
	PolyLineGraph2DHMapImpl plg;
	const GraphAdjacencySetUndirectedNoType &gs = p.first;
	const vector<glm::vec2> &nodes_data = p.second;
	const vector<std::set<ulong>> adjacencies = gs.get_adjacency_sets();
	ulong plg_start,plg_end;
	ulong cur_start,cur_end;
	ulong plg_cur_start, plg_cur_end;
	vector<bool> processed;
	for(ulong i=0;i<nodes_data.size(); i++)
		processed.push_back(false);
	for(ulong i=0;i<nodes_data.size(); i++)
		if(!processed[i]) {
			// Process node i
			const vector<vector<ulong>> polyline_ids_vec = find_polylines(i, adjacencies);

			// Add polylines
			for(const auto polyline_ids : polyline_ids_vec) {
				cur_start = polyline_ids[0];
				plg_cur_start = plg.get_node_id(nodes_data[cur_start]);
				cur_end = polyline_ids[polyline_ids.size()-1];
				plg_cur_end = plg.get_node_id(nodes_data[cur_end]);

				// Set processed, except hubs
				if(!is_hub(adjacencies,cur_start))
					processed[cur_start] = true;
				if(!is_hub(adjacencies,cur_end))
					processed[cur_end] = true;
				for(ulong intermediate_index = 1; intermediate_index < polyline_ids.size()-1; intermediate_index++)
					processed[polyline_ids[intermediate_index]] = true;

				vector<glm::vec2> cur_polyline_coords = get_polyline_coords(polyline_ids, nodes_data);
				// plg checks for duplicates
				plg.add_polyline(cur_polyline_coords);

			}

			processed[i] = true;
		}

	cout << "Generated PolylineGraph from graph: " << plg.get_nodes_amount() << " nodes, " << plg.get_polylines_amount() << " polylines" << endl;

	return plg;
}

PolyLineGraph2DHMapImpl convert_EdgeGraph_to_PolyLineGraph_simplified(const pair<GraphAdjacencySetUndirectedNoType,vector<glm::vec2>> &p) {
	PolyLineGraph2DHMapImpl plg = convert_EdgeGraph_to_PolyLineGraph(p);
	cout << "Polyline non-simplified number of edges: " << plg.get_segments_list().size() << endl;
	plg.simplify();
	return plg;
}

PolyLineGraph2DHMapImpl convertEdgeImagePolyLineGraph_simplified(const Mat &img, const Vec<unsigned char, 3> &edge_color) {
	pair<GraphAdjacencySetUndirectedNoType,vector<glm::vec2>> p = convertEdgeImagePixelToGraph_NoCycles(img,edge_color);
	return convert_EdgeGraph_to_PolyLineGraph_simplified(p);
}

vector<PolyLineGraph2DHMapImpl> convertEdgeImagesPolyLineGraphs_simplified(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color) {
	vector<PolyLineGraph2DHMapImpl> res;

	for(const auto &img : imgs)
		res.push_back(convertEdgeImagePolyLineGraph_simplified(img,edge_color));

	return res;
}

PolyLineGraph2DHMapImpl convertEdgeImagePolyLineGraph(const Mat &img, const Vec<unsigned char, 3> &edge_color) {
	pair<GraphAdjacencySetUndirectedNoType,vector<glm::vec2>> p = convertEdgeImagePixelToGraph_NoCycles(img,edge_color);
	return convert_EdgeGraph_to_PolyLineGraph(p);
}

vector<glm::vec4> convert_graph_to_edge_list_NoCycles(const pair<GraphAdjacencySetUndirectedNoType,vector<glm::vec2>> &p) {
	const GraphAdjacencySetUndirectedNoType &gs = p.first;
	vector<std::set<ulong>> adjacency_list = gs.get_adjacency_sets();
	const vector<glm::vec2> &nodes_data = p.second;
	vector<glm::vec4> res;
	ulong nodes_amount = gs.get_nodes_num();

	for(ulong start_node=0;start_node<nodes_amount;start_node++) {
		const glm::vec2 start_coords = nodes_data[start_node];
		for(const auto end_node : adjacency_list[start_node])
			if(start_node <= end_node) {
				const glm::vec2 end_coords = nodes_data[end_node];
				res.push_back(glm::vec4(start_coords[0],start_coords[1],end_coords[0],end_coords[1]));
			}
	}

	return res;
}

vector<glm::vec4> convertEdgeImagePixelToSegments_NoCycles(const Mat &img, const Vec<unsigned char, 3> &edge_color) {
	pair<GraphAdjacencySetUndirectedNoType,vector<glm::vec2>> p = convertEdgeImagePixelToGraph_NoCycles(img,edge_color);
	return convert_graph_to_edge_list_NoCycles(p);
}

vector<vector<glm::vec4>> convertEdgeImagePixelToSegments_NoCycles(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color) {
	vector<vector<glm::vec4>> res;
	for(const auto &img: imgs)
		res.push_back(convertEdgeImagePixelToSegments_NoCycles(img,edge_color));
	return res;
}

Mat convertEdgeImagePixelToSegmentsImage_NoCycles(const Mat &img, const Vec<unsigned char, 3> &edge_color) {
	vector<glm::vec4> segments = convertEdgeImagePixelToSegments_NoCycles(img,edge_color);
	cout << "Number of segments: " << segments.size() << endl;
	Mat outimg = get_black_image(img);
	draw_segments_on_image(outimg,segments);
	return outimg;
}

vector<Mat> convertEdgeImagesPixelToSegmentsImages_NoCycles(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color) {
	vector<Mat> res;
	for(const auto &img: imgs)
		res.push_back(convertEdgeImagePixelToSegmentsImage_NoCycles(img,edge_color));
	return res;
}

vector<pair<GraphAdjacencySetUndirectedNoType,vector<glm::vec2>>> convertEdgeImagesPixelToGraphs_NoCycles(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color) {
	vector<pair<GraphAdjacencySetUndirectedNoType,vector<glm::vec2>>> res;
	for(const auto &img: imgs)
		res.push_back(convertEdgeImagePixelToGraph_NoCycles(img,edge_color));
	return res;
}

vector<glm::vec4> convertEdgeImagePixelToSegments_PolyLineGraph_simplified(const Mat &img, const Vec<unsigned char, 3> &edge_color) {
	PolyLineGraph2DHMapImpl plg = convertEdgeImagePolyLineGraph_simplified(img,edge_color);
	return plg.get_segments_list();
}

Mat convertEdgeImagePixelToSegmentsImage_PolyLineGraph_simplified(const Mat &img, const Vec<unsigned char, 3> &edge_color) {
	vector<glm::vec4> segments = convertEdgeImagePixelToSegments_PolyLineGraph_simplified(img,edge_color);
	cout << "Simplified PolyLine : Number of segments: " << segments.size() << endl;
	Mat outimg = get_black_image(img);
	draw_segments_on_image(outimg,segments);
	return outimg;
}

vector<Mat> convertEdgeImagesPixelToSegmentsImages_PolyLineGraph_simplified(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color) {
	vector<Mat> res;
		for(const auto &img: imgs)
			res.push_back(convertEdgeImagePixelToSegmentsImage_PolyLineGraph_simplified(img,edge_color));
	return res;
}

Mat convertEdgeImagePixelToSegmentsImage_MultiColor_PolyLineGraph_simplified(const Mat &img, const Vec<unsigned char, 3> &edge_color) {
	vector<glm::vec4> segments = convertEdgeImagePixelToSegments_PolyLineGraph_simplified(img,edge_color);
	cout << "Simplified PolyLine : Number of segments: " << segments.size() << endl;
	Mat outimg = get_black_image(img);
	//draw_segments_on_image(outimg,segments);
	draw_segments_on_image_rnd_colors(outimg,segments);
	return outimg;
}

vector<Mat> convertEdgeImagesPixelToSegmentsImages_MultiColor_PolyLineGraph_simplified(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color) {
	vector<Mat> res;
		for(const auto &img: imgs)
			res.push_back(convertEdgeImagePixelToSegmentsImage_MultiColor_PolyLineGraph_simplified(img,edge_color));
	return res;
}

vector<glm::vec4> convertEdgeImagePixelToSegments_PolyLineGraph(const Mat &img, const Vec<unsigned char, 3> &edge_color) {
	PolyLineGraph2DHMapImpl plg = convertEdgeImagePolyLineGraph(img,edge_color);
	return plg.get_segments_list();
}

Mat convertEdgeImagePixelToSegmentsImage_PolyLineGraph(const Mat &img, const Vec<unsigned char, 3> &edge_color) {
	vector<glm::vec4> segments = convertEdgeImagePixelToSegments_PolyLineGraph(img,edge_color);
	cout << "Simplified PolyLine : Number of segments: " << segments.size() << endl;
	Mat outimg = get_black_image(img);
	draw_segments_on_image(outimg,segments);
	return outimg;
}

vector<Mat> convertEdgeImagesPixelToSegmentsImages_PolyLineGraph(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color) {
	vector<Mat> res;
		for(const auto &img: imgs)
			res.push_back(convertEdgeImagePixelToSegmentsImage_PolyLineGraph(img,edge_color));
	return res;
}

Mat convertEdgeImagePixelToSegmentsImage_MultiColor_PolyLineGraph(const Mat &img, const Vec<unsigned char, 3> &edge_color) {
	vector<glm::vec4> segments = convertEdgeImagePixelToSegments_PolyLineGraph(img,edge_color);
	cout << "Simplified PolyLine : Number of segments: " << segments.size() << endl;
	Mat outimg = get_black_image(img);
	//draw_segments_on_image(outimg,segments);
	draw_segments_on_image_rnd_colors(outimg,segments);
	return outimg;
}

Mat convertEdgeImagePixelToSegmentsImage_MultiColorPolyLines_PolyLineGraph_simplified(const Mat &img, const Vec<unsigned char, 3> &edge_color) {
	PolyLineGraph2DHMapImpl plg = convertEdgeImagePolyLineGraph_simplified(img,edge_color);
	vector<vector<glm::vec4>> vec_segments = plg.get_segments_grouped_by_polyline();
	//vector<glm::vec4> segments = convertEdgeImagePixelToSegments_PolyLineGraph(img,edge_color);
	Mat outimg = get_black_image(img);
	for(const auto &segments : vec_segments)
		draw_segments_on_image_rnd_color(outimg,segments);
	return outimg;
}

vector<Mat> convertEdgeImagesPixelToSegmentsImages_MultiColorPolyLines_PolyLineGraph_simplified(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color) {
	vector<Mat> res;
		for(const auto &img: imgs)
			res.push_back(convertEdgeImagePixelToSegmentsImage_MultiColorPolyLines_PolyLineGraph_simplified(img,edge_color));
	return res;
}

vector<glm::vec2> find_2connections_nodes(const PolyLineGraph2D &plg) {
	vector<glm::vec2> res;

	vector<glm::vec2> node_coords = plg.get_nodes_coords();
	vector<vector<ulong>> connections = plg.get_connections();
	vector<PolyLineGraph2D::polyline> polylines = plg.get_polylines();

	for(ulong node_id = 0; node_id < connections.size(); node_id++)
		if(plg.is_valid_node(node_id))
			if(connections[node_id].size() == 2) {
/*				cout << "Found 2-connect node " << node_id << " with coords " << node_coords[node_id] << " : ";
				print_vector_vec2(polylines[connections[node_id][0]].polyline_coords);
				cout << " / ";
				print_vector_vec2(polylines[connections[node_id][1]].polyline_coords);
				cout << endl;*/
				res.push_back(node_coords[node_id]);
			}

	return res;
}

vector<glm::vec2> find_loops(const PolyLineGraph2D &plg) {
	vector<glm::vec2> res;

	vector<glm::vec2> node_coords = plg.get_nodes_coords();
	vector<PolyLineGraph2D::polyline> polylines = plg.get_polylines();

	for(ulong i=0; i < polylines.size(); i++)
		if(plg.is_valid_polyline(i)) {
			const PolyLineGraph2D::polyline &p = polylines[i];
			if(p.start == p.end) {
/*				cout << "Found loop on node " << p.start << ", with length " << p.polyline_coords.size() << " with coords ";
				print_vector_vec2(p.polyline_coords);
				cout << endl;*/
				res.push_back(node_coords[p.start]);
			}
		}

	return res;
}

Mat convertEdgeImagePixelToSegmentsImage_MultiColorComponents_PolyLineGraph_simplified(const Mat &img, const Vec<unsigned char, 3> &edge_color) {
	PolyLineGraph2DHMapImpl plg = convertEdgeImagePolyLineGraph_simplified(img,edge_color);

	return draw_MultiColorComponents_PolyLineGraph_simplified(img, plg);
}

vector<Mat> convertEdgeImagesPixelToSegmentsImages_MultiColorComponents_PolyLineGraph_simplified(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color) {
	vector<Mat> res;

	for(const auto &img: imgs)
		res.push_back(convertEdgeImagePixelToSegmentsImage_MultiColorComponents_PolyLineGraph_simplified(img,edge_color));

	return res;
}

vector<vector<glm::vec4>> convertEdgeImagePixelToSegments_PolyLineGraph_simplified(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color) {
	vector<PolyLineGraph2DHMapImpl> plgs = convertEdgeImagesPolyLineGraphs_simplified(imgs,edge_color);
	vector<vector<glm::vec4>> vec_segments;
	for(const auto &plg : plgs)
		vec_segments.push_back(plg.get_segments_list());
	return vec_segments;
}

Mat convertEdgeImagePixelToSegmentsImage_MultiColorComponents_PolyLineGraph_optimized(const Mat &img, const Vec<unsigned char, 3> &edge_color) {
	PolyLineGraph2DHMapImpl plg = convertEdgeImagePolyLineGraph_optimized(img,edge_color);

	return draw_MultiColorComponents_PolyLineGraph_simplified(img, plg);
}

vector<Mat> convertEdgeImagesPixelToSegmentsImages_MultiColorComponents_PolyLineGraph_optimized(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color) {
	vector<Mat> res;

	for(const auto &img: imgs)
		res.push_back(convertEdgeImagePixelToSegmentsImage_MultiColorComponents_PolyLineGraph_optimized(img,edge_color));

	return res;
}

PolyLineGraph2DHMapImpl convert_EdgeGraph_to_PolyLineGraph_optimized(const pair<GraphAdjacencySetUndirectedNoType,vector<glm::vec2>> &p) {
	PolyLineGraph2DHMapImpl plg = convert_EdgeGraph_to_PolyLineGraph(p);
	cout << "Polyline non-simplified number of edges: " << plg.get_segments_list().size() << endl;
	plg.optimize();
	return plg;
}

/**
 * Adds:
 * - polyline simplification
 * - direct extreme connect
 */
PolyLineGraph2DHMapImpl convertEdgeImagePolyLineGraph_optimized(const Mat &img, const Vec<unsigned char, 3> &edge_color) {
	pair<GraphAdjacencySetUndirectedNoType,vector<glm::vec2>> p = convertEdgeImagePixelToGraph_NoCycles(img,edge_color);
	return convert_EdgeGraph_to_PolyLineGraph_optimized(p);
}

vector<PolyLineGraph2DHMapImpl> convert_edge_images_to_optimized_polyline_graphs(const vector<Mat> &imgs, const Vec<unsigned char, 3> &edge_color) {
	vector<PolyLineGraph2DHMapImpl> res;

	for(const auto &img : imgs)
		res.push_back(convertEdgeImagePolyLineGraph_optimized(img,edge_color));

	return res;
}
