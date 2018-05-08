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


#ifndef INCLUDE_EDGEGRAPH3D_PLGS_POLYLINE_GRAPH_3D_HPP_
#define INCLUDE_EDGEGRAPH3D_PLGS_POLYLINE_GRAPH_3D_HPP_

#include <utility>
#include <vector>

#include "datatypes.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"
#include "glm.hpp"


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

#define TOP_FILTER_BY_POLYLINESMOOTHLENGTH 0.95

#define MAXIMUM_LINEARIZABILITY_DISTANCE 0.01
#define MAXIMUM_LINEARIZABILITY_DISTANCE_SQ (MAXIMUM_LINEARIZABILITY_DISTANCE*MAXIMUM_LINEARIZABILITY_DISTANCE)

class PolyLineGraph3D {
public:

	struct polyline {
		ulong start;
		ulong end;
		vector<glm::vec3> polyline_coords;
		float length;

		polyline();
		polyline(const ulong s, const ulong e, const vector<glm::vec3> &pcs);
		polyline(const ulong &s, const ulong &e, const vector<glm::vec3> &pcs, float length);
		/**
		 * return pair(begin,end)
		 * 	->	return iterators if begin_id is start, reverse iterators otherwise
		 */
		//pair<vector<glm::vec3>::iterator,vector<glm::vec3>::iterator> get_iterator(ulong begin_id);
		void simplify();
		void simplify(const float max_linearizability_dist);

		void update_length();
		//vector<glm::vec6> get_segments_list();
		vector<vec6> get_segments_list();

		bool polyline::operator==(const polyline& a) const;

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

		bool is_start(const ulong node_id);
		bool is_end(const ulong node_id);

		bool connects(const ulong a);
		bool connects(const ulong a, const ulong b);

		void fragment(const float maxlen);

		float get_maxlength();

		struct pl_point {
			const ulong segment_index; // 5 means that the match is found on the segment that goes from node 5 to node 6 in the polyline coords
			const glm::vec3 coords;
			pl_point(const ulong segment_index, const glm::vec3 &coords);
		};

	};

	void fragment(const float maxlen);

	vector<polyline> polylines;

	// for each node, vector of ids of connected polylines
	vector<vector<ulong>> connections;

	vector<bool> visited_nodes;
	ulong real_nodes_amount; // No removed nodes
	ulong nodes_amount;
	ulong next_node_id;
	vector<glm::vec3> nodes_coords;

	vector<pair<vector<glm::vec2>,vector<int>>> observations;

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

	struct plg_point {
		ulong polyline_id;
		ulong segment_index; // 5 means that the match is found on the segment that goes from node 5 to node 6 in the polyline coords
		glm::vec3 coords;

		plg_point();
		plg_point(const ulong input_polyline_id, const ulong input_segment_index, const glm::vec3 input_coords);
		plg_point(const pair<ulong,ulong> &input_polyline_segment_index, const glm::vec3 input_coords);
	};

	PolyLineGraph3D();



	virtual ulong get_node_id(const glm::vec3 &p_coords) = 0;
	virtual void add_polyline(const vector<glm::vec3> &polyline_coords) = 0;
	vector<vec6> get_segments_list();
	vector<vector<vec6>> get_segments_list_by_polyline();
	vector<vec6> get_segments_list(set3dpoints include_only);
	glm::vec3 get_node_coords(ulong node_id);
	vector<glm::vec3> get_nodes_coords();
	vector<vector<ulong>> get_connections();
	vector<polyline> get_polylines();
	ulong get_nodes_amount();
	ulong get_polylines_amount();
	void simplify();
	void simplify(const float max_linearizability_dist);
	void optimize();
	virtual void add_direct_connection(const ulong start, const ulong end) = 0;
	//virtual void connect_close_extremes() = 0;
	ulong get_real_nodes_amount();
	bool is_valid_node(const ulong node_id);
	bool is_valid_polyline(const ulong polyline_id);
	void remove_polyline(ulong polyline_id);
	bool has_loop(ulong node_id);
	void set_observations(ulong node_id, vector<glm::vec2> projections, vector<int> cam_ids);

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
	vector<ulong> PolyLineGraph3D::get_nodes_with_loops();

	vector<glm::vec3> get_hub_nodes_coords();
	vector<glm::vec3> get_extreme_nodes_coords();
	vector<glm::vec3> get_loopnodes_coords();
	vector<glm::vec3> get_nodes_with_loops_coords();

	bool is_connected_node_inside_radius(const ulong start_node,const ulong end_node, const float detection_radius);
	bool is_connected_polyline(ulong polyline_id_start,ulong polyline_id_end);
	bool is_connected_polyline(ulong polyline_id_start,ulong polyline_id_end, ulong max_jumps);
	bool is_connected_polyline_inside_radius(ulong polyline_id_start, ulong polyline_id_end, const float detection_radius);

	void print_stats();

protected:
	void add_node_coords(const glm::vec3 &p_coords);
	virtual ~PolyLineGraph3D() {};
	ulong get_next_node_id();
	void increment_real_nodes_amount();
	void decrement_real_nodes_amount();
	void invalidate_node(ulong node_id);
	void remove_connection(const ulong node_id, const ulong polyline_id);
	bool is_connected_node(ulong start_node,ulong end_node);
	bool is_connected_node(ulong start_node,ulong end_node, ulong max_jumps);

	vector<ulong> valid_nodes_index(ulong &valid_nodes_amount);
	ulong get_amount_valid_polylines();

	PolyLineGraph3D(const vector<polyline> &polylines, const vector<vector<ulong>> &connections, const vector<bool> &visited_nodes, const ulong &real_nodes_amount, const ulong &nodes_amount, const ulong &next_node_id, const vector<glm::vec3> &nodes_coords);

private:

};

#endif /* INCLUDE_EDGE_MATCHER_EDGES_OUTPUT_CONVERTER_POLYLINEGRAPH2D_HPP_ */
