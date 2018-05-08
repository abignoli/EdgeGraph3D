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


#ifndef INCLUDE_EDGEGRAPH3D_PLGS_POLYLINE_GRAPH_3D_HMAP_IMPL_HPP_
#define INCLUDE_EDGEGRAPH3D_PLGS_POLYLINE_GRAPH_3D_HMAP_IMPL_HPP_

#include <boost/serialization/base_object.hpp>
#include <stddef.h>
#include <string>
#include <utility>
#include <vector>
#include <unordered_map>
#include <fstream>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "polyline_graph_3d.hpp"
#include "datatypes.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"
#include "serialization_utilities.hpp"
#include "glm.hpp"


namespace boost {
namespace serialization {
class access;
} /* namespace serialization */
} /* namespace boost */

using namespace std;

struct KeyFuncs3d
{
    size_t operator()(const glm::vec3& k)const
    {
        return std::hash<int>()(k.x) ^ std::hash<int>()(k.y) ^ std::hash<int>()(k.z);
    }

    bool operator()(const glm::vec3& a, const glm::vec3& b)const
    {
            return a.x == b.x && a.y == b.y && a.z == b.z;
    }
};

typedef unordered_map<glm::vec3,ulong,KeyFuncs3d,KeyFuncs3d> pointmap3dtype;

class PolyLineGraph3DHMapImpl: public PolyLineGraph3D {
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        // serialize base class information
        ar & boost::serialization::base_object<PolyLineGraph3D>(*this);
        ar & point_map;
    }


	bool is_duplicate(const polyline &pl);
	void internal_add_polyline(const polyline &pl);
	void invalidate_node(ulong node_id);
	void remove_invalid_polylines();
public:
	PolyLineGraph3DHMapImpl();
	ulong get_node_id(const glm::vec3 &p_coords);
	void remap_node(const ulong node_id, const glm::vec3 &new_coords);
	void remap_node(const glm::vec3 &old_coords, const glm::vec3 &new_coords);
	void add_polyline(const vector<glm::vec3> &polyline_coords);
	void add_direct_connection(const ulong start, const ulong end);
	void add_direct_connection(const glm::vec3 &start_coords, const glm::vec3 &end_coords);
	void add_direct_connection(const glm::vec3 &start_coords, const glm::vec3 &end_coords, pair<ulong, ulong> &out_ids);
	void filter_nodes(const set3dpoints &inliers);
	void remove_polylines_with_longsegments(const float toplength_ratio);

	pointmap3dtype point_map;
	PolyLineGraph3DHMapImpl(const vector<polyline> &polylines, const vector<vector<ulong>> &connections, const vector<bool> &visited_nodes, const ulong &real_nodes_amount, const ulong &nodes_amount, const ulong &next_node_id, const vector<glm::vec3> &nodes_coords, const pointmap3dtype &point_map);

	~PolyLineGraph3DHMapImpl();
};

void serialize_plg(const PolyLineGraph3DHMapImpl &plg, const string plg_path);

PolyLineGraph3DHMapImpl deserialize_3dplg(const string plg_path);

#endif /* INCLUDE_EDGE_MATCHER_EDGES_OUTPUT_CONVERTER_POLYLINEGRAPH2DHMAPIMPL_HPP_ */
