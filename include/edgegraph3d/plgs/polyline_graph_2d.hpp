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


#ifndef INCLUDE_EDGEGRAPH3D_PLGS_POLYLINE_GRAPH_2D_HPP_
#define INCLUDE_EDGEGRAPH3D_PLGS_POLYLINE_GRAPH_2D_HPP_

#include <iostream>
#include <set>
#include <utility>
#include <vector>

#include "datatypes.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"
#include "glm.hpp"


struct ray2d;

namespace boost {
namespace serialization {
class access;
} /* namespace serialization */
} /* namespace boost */

using namespace std;

#define DIRECT_CONNECTION_EXTREMES_MAXDIST 6

#define DIRECT_CONNECTION_EXTREMES_FOLLOWING_DIRECTION_MAXDIST 15
#define DIRECT_CONNECTION_EXTREMES_FOLLOWING_DIRECTION_MINCOS 0.707

#define INVALID_POINT_COORDS -1
#define INVALID_POLYLINE_LENGTH -1

#define SMOOTHSECTION_MAX_ANGLE 45
#define SMOOTHSECTION_MIN_COS 0.707

#define TOP_FILTER_BY_POLYLINESMOOTHLENGTH 0.82

#define MAXIMUM_LINEARIZABILITY_DISTANCE 1.0
#define MAXIMUM_LINEARIZABILITY_DISTANCE_SQ (MAXIMUM_LINEARIZABILITY_DISTANCE*MAXIMUM_LINEARIZABILITY_DISTANCE)

#define POLYLINE_NEXT_BY_LINE_INTERSECTION_MAX_QUASIPARALLEL_ANGLE 15
#define POLYLINE_NEXT_BY_LINE_INTERSECTION_MAX_QUASIPARALLEL_ANGLE_COS 0.965
#define POLYLINE_NEXT_BY_LINE_INTERSECTION_MAX_QUASIPARALLEL_DIST 5

#define POLYLINE_RAY_INTERSECT_APPROX_MAX_DIST 1.0
#define POLYLINE_RAY_INTERSECT_APPROX_MAX_DISTSQ (POLYLINE_RAY_INTERSECT_APPROX_MAX_DIST*POLYLINE_RAY_INTERSECT_APPROX_MAX_DIST)

#define PROLONG_EXTREME_MIN_SEGMENT_LENGTH 5
#define PROLONG_EXTREME_MAX_PROLONGATION 5

class PolyLineGraph2D {
public:

	struct polyline {

		struct pl_point {
			ulong segment_index; // 5 means that the match is found on the segment that goes from node 5 to node 6 in the polyline coords
			glm::vec2 coords;
			template<class Archive>
			void serialize(Archive & ar, const unsigned int version)
			{
			    ar & segment_index;
			    ar & coords;
			}
			pl_point();
			pl_point(const ulong segment_index, const glm::vec2 &coords);
			bool operator==(const pl_point& a) const;
			//pl_point(const pl_point &plp);
			//PolyLineGraph2D::polyline::pl_point& PolyLineGraph2D::polyline::pl_point::operator=(const PolyLineGraph2D::polyline::pl_point&)
		};

		struct pl_interval {
			pl_point start;
			pl_point end;
			template<class Archive>
			void serialize(Archive & ar, const unsigned int version)
			{
			    ar & start;
			    ar & end;
			}
			pl_interval();
			pl_interval(const pl_point start, const pl_point end);
		};
		bool interval_contains_plp(const pl_interval &pli, const pl_point &plp);

		ulong start;
		ulong end;
		vector<glm::vec2> polyline_coords;
		float length;

		polyline();
		polyline(const ulong s, const ulong e, const vector<glm::vec2> &pcs);
		polyline(const ulong &s, const ulong &e, const vector<glm::vec2> &pcs, float length);
		~polyline();

		vector<pair<ulong,ulong>> get_intersectedcells_2dmap_vec(const float cell_dim);
		set<pair<ulong,ulong>> get_intersectedcells_2dmap_set(const float cell_dim);

		/**
		 * return pair(begin,end)
		 * 	->	return iterators if begin_id is start, reverse iterators otherwise
		 */
		//pair<vector<glm::vec2>::iterator,vector<glm::vec2>::iterator> get_iterator(ulong begin_id);
		void simplify();
		void simplify(const float max_linearizability_dist);

		void update_length();
		vector<glm::vec4> get_segments_list();

		bool operator==(const polyline& a) const;

		static polyline merge_polylines(const polyline& p1,const polyline& p2);

		void clear_coords();
		void invalidate();

		bool is_loop();

		// Given start returns end and viceversa
		ulong get_other_end(ulong extreme);

		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
		    ar & start;
		    ar & end;
		    ar & polyline_coords;
		    ar & length;
		}

		ulong get_extreme_id(const pl_point &plp);
		bool is_start(const ulong node_id);
		bool is_start(const pl_point &plp);
		bool is_end(const ulong node_id);
		bool is_end(const pl_point &plp);

		glm::vec2 get_extreme_coordinates(const ulong extreme_id);
		pl_point get_start_plp();
		pl_point get_end_plp();
		pl_point get_extreme_plp(const ulong extreme_id);
		pl_point get_extreme_plp(const ulong extreme_id, bool &valid);

		bool connects(const ulong a);
		bool connects(const ulong a, const ulong b);

		float compute_distancesq(const glm::vec2 &p_coords, ulong &closest_segm, glm::vec2 &projection); // Compute distance from point
		float compute_distancesq(const glm::vec2 &p_coords, pl_point &plp); // Compute distance from point

		/**
		 * Returns the length of the longest polyline portion considered smooth
		 *
		 * A polyline portion is considered smooth when each segment's direction
		 * produces an angle no bigger than MAX_ANGLE (i.e. cos > MIN_COS).
		 */
		float compute_max_smooth_length();
		// get_iterator_from_ppline

		bool has_point(const pl_point &plp);

		// Input: start or end id
		// Output: direction and length of last segment
		pair<glm::vec2,float> get_extreme_direction_length(const ulong extreme);

		// Input: start or end id
		// Output: direction and length of last portion of polyline with given lengthsq
		pair<glm::vec2,float> get_extreme_direction_length_given_length(const ulong extreme, const float length);

		pair<vector<glm::vec2>,vector<glm::vec2>> split(const pl_point &plp);

		vector<pl_point> intersect_segment(const glm::vec4 &segment);
		vector<pl_point> intersect_line(const glm::vec3 &line);
		vector<pl_point> intersect_ray(const ray2d &r);

		void first_intersect_ray(const ray2d &r, pl_point &intersection, float &distance, bool &found);
		// Approximate intersection to closest polyline coords if close enough (POLYLINE_RAY_INTERSECT_APPROX_MAX_DIST)
		void first_intersect_ray_approx(const ray2d &r, pl_point &intersection, float &distance, bool &found);

		// direction can be either start or end
		//void next_pl_point_radius(const ulong init_segment_index, const glm::vec2 init_coords, const ulong direction, ulong &next_segment_index, glm::vec2 &next_coords);
		pl_point next_pl_point_by_distance(const PolyLineGraph2D::polyline::pl_point init_plp, const ulong direction, const float distance, bool &reached_polyline_extreme);

		// direction can be either start or end
		//void next_pl_point_radius(const ulong init_segment_index, const glm::vec2 init_coords, const ulong direction, ulong &next_segment_index, glm::vec2 &next_coords);
		pl_point next_pl_point_by_length(const PolyLineGraph2D::polyline::pl_point init_plp, const ulong direction, const float length, bool &reached_polyline_extreme);

		// Get next point by distance from specified coords, starting from specified polyline extreme
		pl_point next_pl_point_by_distance(const ulong polyline_extreme, const glm::vec2 coords, const float distance, bool &reached_polyline_extreme);

		vector<pl_point> next_pl_points_by_distance(const pl_point init_plp, const ulong direction, const float distance);
		vector<pl_point> split_equal_size_intervals(const ulong starting_extreme, const float distance);

		//pl_point next_pl_point_by_length(const pl_point init_plp, const ulong direction, const float length);

		void next_pl_point_by_line_intersection(const pl_point init_plp, const ulong direction, const glm::vec3 &line, PolyLineGraph2D::polyline::pl_point &next, bool &found_quasiparallel_segment, PolyLineGraph2D::polyline::pl_point &next_before_quasiparallel_segment, bool &reached_polyline_extreme, bool &found);
		//void void PolyLineGraph2D::polyline::closest_next_pl_point_by_line_intersection_nodirection(const pl_point init_plp, const glm::vec3 &line, PolyLineGraph2D::polyline::pl_point &next, bool &found_quasiparallel_segment, bool &found);
		void next_pl_point_by_line_intersection_bounded_distance(const pl_point init_plp, const ulong direction, const glm::vec3 &line, const float min_dist, const float max_dist, PolyLineGraph2D::polyline::pl_point &next, bool &found_quasiparallel_segment, PolyLineGraph2D::polyline::pl_point &next_before_quasiparallel_segment, bool &reached_polyline_extreme, bool &bounded_distance_violated, bool &found);
		ulong get_amount_of_segments();

		//void find_3dpoint_ppline(const new_3dpoint_plgp_matches &p3d, bool &valid);
	};

	vector<polyline> polylines;

	// for each node, vector of ids of connected polylines
	vector<vector<ulong>> connections;

	vector<bool> visited_nodes;
	ulong real_nodes_amount; // No removed nodes
	ulong nodes_amount;
	ulong next_node_id;
	vector<glm::vec2> nodes_coords;

    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & polylines;
        ar & connections;
        ar & visited_nodes;
        ar & real_nodes_amount;
        ar & nodes_amount;
        ar & next_node_id;
        ar & nodes_coords;
    }

/*	struct polyline_store {
		vector<polyline> polylines;
		stack<ulong> invalidated_ids;
		// Add a polyline and return its ID
		ulong add_polyline(polyline &p);
		void remove_polyline(const ulong polyline_id);
	};

	struct node_store {
		vector<vector<ulong>> connections;
		vector<glm::vec2> nodes_coords;
		stack<ulong> invalidated_ids;
		// Add a node and return its ID
		ulong add_node(const glm::vec2 &coords);
		void remove_node(const ulong node_id);
		void remove_node_polyline_connection(const ulong node_id, const ulong polyline_id);
	};*/

	struct plg_point {
		ulong polyline_id;
/*		ulong segment_index; // 5 means that the match is found on the segment that goes from node 5 to node 6 in the polyline coords
		glm::vec2 coords;*/
		polyline::pl_point plp;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
		    ar & polyline_id;
		    ar & plp;
		}
		plg_point();
		plg_point(const ulong input_polyline_id, const ulong input_segment_index, const glm::vec2 input_coords);
		plg_point(const ulong input_polyline_id, const polyline::pl_point &plp);
		plg_point(const pair<ulong,ulong> &input_polyline_segment_index, const glm::vec2 input_coords);
		bool operator==(const plg_point& a) const;
	};

	struct plgs_point {
		ulong plg_id;
		plg_point plgp;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
		    ar & plg_id;
		    ar & plgp;
		}
		plgs_point();
		plgs_point(const ulong plg_id, const plg_point &plgp);
		bool operator==(const plgs_point& a) const;

	};

	struct plg_node_components {
		vector<ulong> nodes_component_id;
		vector<set<ulong>> components_nodes_ids;
		plg_node_components(vector<ulong> &nodes_component_id,vector<set<ulong>> &components_nodes_ids);
		~plg_node_components();
	};

	struct plg_polyline_components {
		vector<ulong> polylines_component_id;
		vector<set<ulong>> components_polylines_ids;
		plg_polyline_components(vector<ulong> &polylines_component_id,vector<set<ulong>> &components_polylines_ids);
		~plg_polyline_components();
	};

	struct plg_components {
		plg_node_components nc;
		plg_polyline_components pc;
		plg_components(plg_node_components &nc, plg_polyline_components &pc);
		plg_components(vector<ulong> &nodes_component_id,vector<set<ulong>> &components_nodes_ids, vector<ulong> &polylines_component_id,vector<set<ulong>> &components_polylines_ids);
		~plg_components();
	};

	PolyLineGraph2D();

	virtual ulong get_node_id(const glm::vec2 &p_coords) = 0;
	virtual void add_polyline(const vector<glm::vec2> &polyline_coords) = 0;
	vector<glm::vec4> get_segments_list();
	vector<vector<glm::vec4>> get_segments_grouped_by_polyline();
	vector<pair<ulong,vector<glm::vec4>>> get_segments_grouped_by_polyline_with_polyline_ids();
	vector<vector<glm::vec4>> get_segments_grouped_by_component();
	glm::vec2 get_node_coords(ulong node_id);
	vector<glm::vec2> get_nodes_coords();
	vector<vector<ulong>> get_connections();
	vector<polyline> get_polylines();
	ulong get_nodes_amount();
	ulong get_polylines_amount();
	void simplify();
	void simplify(const float max_linearizability_dist);
	void optimize();
	virtual void add_direct_connection(const ulong start, const ulong end) = 0;
	virtual void connect_close_extremes() = 0;
	ulong get_real_nodes_amount();
	bool is_valid_node(const ulong node_id);
	bool is_valid_polyline(const ulong polyline_id);
	void remove_polyline(ulong polyline_id);
	bool has_loop(ulong node_id);
	vector<ulong> compute_components_node_ids();
	pair<vector<ulong>,vector<set<ulong>>> compute_components(); // nodes
	plg_components compute_plg_components();
	pair<pair<vector<ulong>,vector<set<ulong>>>,pair<vector<ulong>,vector<set<ulong>>>> compute_components_with_polylines();

	void filter_components_by_polylinesmoothlength();

	vector<pair<ulong,vector<polyline::pl_point>>> intersect_polylines(const glm::vec4 &segment);

	/* Hub, either:
	 * - >= 3 connected polylines
	 * - 2 connected polylines, provided 1 polyline is a loop
	 */
	bool is_hub(ulong node_id);

	/* Extreme, either:
	 * - 1 connected polyline (not a loop)
	 */
	bool is_extreme(ulong node_id);

	/* LoopNode, either:
	 * - 1 connected polyline (a loop)
	 */
	bool is_loopnode(ulong node_id);

	vector<ulong> get_hub_nodes();
	vector<ulong> get_extreme_nodes();
	vector<ulong> get_loopnodes();
	vector<ulong> get_nodes_with_loops();

	pair<vector<ulong>,vector<glm::vec2>> get_extreme_nodes_ids_and_coords();
	tuple<vector<ulong>,vector<glm::vec2>,vector<glm::vec2>> get_extreme_nodes_ids_and_coords_and_direction();

	vector<glm::vec2> get_hub_nodes_coords();
	vector<glm::vec2> get_extreme_nodes_coords();
	vector<glm::vec2> get_loopnodes_coords();
	vector<glm::vec2> get_nodes_with_loops_coords();

	bool is_connected_node_inside_radius(const ulong start_node,const ulong end_node, const float detection_radius);
	bool is_connected_polyline(ulong polyline_id_start,ulong polyline_id_end);
	bool is_connected_polyline(ulong polyline_id_start,ulong polyline_id_end, ulong max_jumps);
	bool is_connected_polyline_inside_radius(ulong polyline_id_start, ulong polyline_id_end, const float detection_radius);

	/*


		 * CPF - Closest Polyline Finder
		 *
		 * The following code has the function of enabling quick determination
		 * of the closest polyline/segments to a given point
		 *
		 * point coordinates -> closest_polyline


		bool cpf_init = false;
		const float max_detection_dist = DETECTION_CORRESPONDENCES_RADIUS;
		vector<pair<ulong,ulong>> **polymap;
		void init_cpf(const ulong img_rows, const ulong img_cols);
		vector<tuple<float,ulong,ulong,glm::vec2>> cpf_find_nearby_polylines(const glm::vec2 &coords); // Find polylines within DETECTION_CORRESPONDENCES_RADIUS from given coords
	*/

	void intersect_ray_first_polyline_within_dist(const ray2d &r, const float max_distance, plg_point &plgp, float &distance, bool &found);
	void intersect_ray_first_polyline_within_dist(const ray2d &r, const ulong polyline_to_exclude, const float max_distance, plg_point &plgp, float &distance, bool &found);

	float cpf_find_unbound(const glm::vec2 &coords, plg_point &plgp);
	float cpf_find_unbound(const glm::vec2 &coords, ulong &closest_polyline, ulong &closest_segment_on_polyline, glm::vec2 &projection_coords);

	// Output: for each close polyline : <polyline_id,segment_index,closest_point_coords,distance>
	vector<tuple<ulong,ulong,glm::vec2,float>> cpf_find_within_radius(const glm::vec2 &coords, const float max_dist);
	//vector<tuple<ulong,ulong,glm::vec2,float>> cpf_find_within_radius_component(const glm::vec2 &coords, const float max_dist);
	//vector<tuple<ulong,ulong,glm::vec2,float>> cpf_find_within_radius_parallel(const glm::vec2 &coords, const float max_dist);

/*
	// Output: for each close component : <polyline_id,segment_index,closest_point_coords,distance>
	vector<tuple<ulong,ulong,glm::vec2,float>> cpf_find_within_radius(const glm::vec2 &coords, const float max_dist);
*/

	//void add_polyline(const polyline &pl);
protected:
	void add_node_coords(const glm::vec2 &p_coords);
	ulong get_next_node_id();
	void increment_real_nodes_amount();
	void decrement_real_nodes_amount();
	void invalidate_node(ulong node_id);
	void remove_connection(const ulong node_id, const ulong polyline_id);
	bool is_connected_node(ulong start_node,ulong end_node);
	bool is_connected_node(ulong start_node,ulong end_node, ulong max_jumps);

	PolyLineGraph2D(const vector<polyline> &polylines, const vector<vector<ulong>> &connections, const vector<bool> &visited_nodes, const ulong &real_nodes_amount, const ulong &nodes_amount, const ulong &next_node_id, const vector<glm::vec2> &nodes_coords);
	virtual ~PolyLineGraph2D();
private:

};

typedef std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> new_3dpoint_plgp_matches;
void update_new_3dpoint_plgp_matches(new_3dpoint_plgp_matches &pt, const int new_view, const PolyLineGraph2D::plg_point &new_observation, const glm::vec3 &new_coords);

/**
 * new_3dpoint_and_sides_plgp_matches
 * tuple<pair<vector<new_3dpoint_plgp_matches>,vector<ulong>>,new_3dpoint_plgp_matches,pair<vector<new_3dpoint_plgp_matches>,vector<ulong>>>
 * 			direction1: (3d points,				directions) , central 3d point			,direction2: (3d points,				directions)
 */
typedef tuple<pair<vector<new_3dpoint_plgp_matches>,vector<ulong>>,new_3dpoint_plgp_matches,pair<vector<new_3dpoint_plgp_matches>,vector<ulong>>> new_3dpoint_and_sides_plgp_matches;

std::ostream &operator<< (std::ostream &out, const PolyLineGraph2D::plg_point &plgp);

vector<new_3dpoint_plgp_matches> new_3dpoint_and_sides_plgp_matches_to_vector(const new_3dpoint_and_sides_plgp_matches &p3d_with_sides);

new_3dpoint_and_sides_plgp_matches make_new_3dpoint_and_sides_plgp_matches(vector<new_3dpoint_plgp_matches> valid_points_direction1,vector<ulong> directions1,new_3dpoint_plgp_matches central_point,vector<new_3dpoint_plgp_matches> valid_points_direction2,vector<ulong> directions2);

vector<PolyLineGraph2D::plg_point> convert_vec_pl_point_to_plg_point(const vector<PolyLineGraph2D::polyline::pl_point> &vpl, const ulong polyline_id);

vector<glm::vec2> convert_vec_pl_point_to_vec2(const vector<PolyLineGraph2D::polyline::pl_point> &vpl);

//vector<pair<glm::vec2,glm::vec2>> find_closest_pairs_with_max_dist(const vector<glm::vec2> &p2ds, const float max_dist);
vector<pair<ulong,ulong>> find_closest_pairs_with_max_dist(pair<vector<ulong>,vector<glm::vec2>> nodes_data, const float max_dist);
vector<pair<ulong,ulong>> find_closest_pairs_with_max_dist_following_direction(tuple<vector<ulong>,vector<glm::vec2>,vector<glm::vec2>> nodes_data, const float max_dist, const float min_cos);

bool one_or_zero_correspondences(const pair<PolyLineGraph2D::plg_point,std::vector<std::vector<PolyLineGraph2D::plg_point>>> &p);

int amount_of_total_2d_correspondences(const pair<PolyLineGraph2D::plg_point,std::vector<std::vector<PolyLineGraph2D::plg_point>>> &p);

vector<glm::vec2> convert_vecplgp_to_vec2(const std::vector<PolyLineGraph2D::plg_point> &plgps);

pair<vector<glm::vec2>,vector<int>> convert_plgpoint_correspondences(const vector<int> &cams, const std::vector<std::vector<PolyLineGraph2D::plg_point>> &potential_match_correspondences);

pair<vector<PolyLineGraph2D::plg_point>,vector<int>> convert_plgpoint_correspondences_plgp(const vector<int> &cams, const std::vector<std::vector<PolyLineGraph2D::plg_point>> &potential_match_correspondences);

#endif /* INCLUDE_EDGEGRAPH3D_PLGS_POLYLINE_GRAPH_2D_HPP_ */
