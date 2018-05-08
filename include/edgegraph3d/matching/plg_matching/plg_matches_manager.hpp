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


#ifndef INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_PLG_MATCHES_MANAGER_HPP_
#define INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_PLG_MATCHES_MANAGER_HPP_

#include <omp.h>
#include <stddef.h>
#include <set>
#include <string>
#include <vector>

#include "glm.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/set.hpp>

#include "polyline_graph_2d.hpp"
#include "polyline_graph_2d_hmap_impl.hpp"
#include "polyline_graph_3d_hmap_impl.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"
#include "global_switches.hpp"
#include "serialization_utilities.hpp"
#include "tuple_serialization.hpp"

namespace boost {
namespace serialization {
class access;
} /* namespace serialization */
} /* namespace boost */

struct interval_compare {
    bool operator() (const PolyLineGraph2D::polyline::pl_interval& a, const PolyLineGraph2D::polyline::pl_interval& b) const{
    	return a.start.segment_index < b.start.segment_index;
    }
};
typedef set<PolyLineGraph2D::polyline::pl_interval, interval_compare> sorted_pl_intervals_set;

struct KeyFuncs_plgsp_to_p3d_id_map
{
    size_t operator()(const PolyLineGraph2D::plgs_point& k)const
    {
    	return std::hash<int>()(k.plg_id) ^ std::hash<ulong>()(k.plgp.polyline_id) ^ std::hash<ulong>()(k.plgp.plp.segment_index) ^ std::hash<int>()(k.plgp.plp.coords.x) ^ std::hash<int>()(k.plgp.plp.coords.y);
    }

    bool operator()(const PolyLineGraph2D::plgs_point& a, const PolyLineGraph2D::plgs_point& b)const
    {
            return a==b;
    }
};

typedef unordered_map<PolyLineGraph2D::plgs_point,ulong,KeyFuncs_plgsp_to_p3d_id_map,KeyFuncs_plgsp_to_p3d_id_map> plgsp_to_p3d_id_map;

struct KeyFuncs3dtoplgps
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

typedef unordered_map<glm::vec3,new_3dpoint_plgp_matches,KeyFuncs3dtoplgps,KeyFuncs3dtoplgps> pointmap3dtoplgpstype;

class PLGMatchesManager {
private:
	omp_lock_t writelock;

	// For each PLG, For each polyline, For each segment: true if there is already a match
	//vector<vector<vector<bool>>> matched_segments;
public:
	friend class boost::serialization::access;
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        // serialize base class information
        ar & plgs;
        ar & plg3d;
        ar & point_matches;
        ar & point_matches_map;
        ar & p3d_to_matches_map;
        ar & matched_polyline_intervals;
    }

	vector<PolyLineGraph2DHMapImpl> &plgs;
	PolyLineGraph3DHMapImpl &plg3d;
	// For each 3D point in plg3d: vector of 2D plgps that generated it
	vector<PolyLineGraph2D::plgs_point> point_matches;
	plgsp_to_p3d_id_map point_matches_map; // maps each plgsp to the id of the 3D point in plg3d
	pointmap3dtoplgpstype p3d_to_matches_map;
	// For each PLG, For each polyline, set of pair(pl_point,pl_point) i.e. start and end of already matched interval
	vector<vector<sorted_pl_intervals_set>> matched_polyline_intervals;


	void add_matched_2dsegment(const int plg_id, const ulong pl_id, const PolyLineGraph2D::polyline::pl_point &pl_a, const PolyLineGraph2D::polyline::pl_point &pl_b);
	void add_matched_3dsegment(const new_3dpoint_plgp_matches &p1, const new_3dpoint_plgp_matches &p2);
	void add_matched_3dpolyline(const vector<new_3dpoint_plgp_matches> &pl);

	bool is_matched(const int plg_id,const ulong polyline_id, const PolyLineGraph2D::polyline::pl_point &plp, PolyLineGraph2D::polyline::pl_interval &pli_containing_point);
	bool is_matched(const int plg_id, const PolyLineGraph2D::plg_point &plgp, PolyLineGraph2D::polyline::pl_interval &pli_containing_point);
	bool is_matched(const PolyLineGraph2D::plgs_point &plgsp, PolyLineGraph2D::polyline::pl_interval &pli_containing_point);

	PolyLineGraph3DHMapImpl& get_plg3d();

	PLGMatchesManager();
	PLGMatchesManager(vector<PolyLineGraph2DHMapImpl> &plgs, PolyLineGraph3DHMapImpl &plg3d);
	PLGMatchesManager(vector<PolyLineGraph2DHMapImpl> &plgs, PolyLineGraph3DHMapImpl &plg3d, vector<PolyLineGraph2D::plgs_point> &point_matches, plgsp_to_p3d_id_map &point_matches_map, pointmap3dtoplgpstype &p3d_to_matches_map, vector<vector<sorted_pl_intervals_set>> &matched_polyline_intervals);
	~PLGMatchesManager();
};

void serialize_plgmm(const PLGMatchesManager &plgmm, const string plgmm_path);

PLGMatchesManager deserialize_plgmm(const string plgmm_path);

#endif /* INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_PLG_MATCHES_MANAGER_HPP_ */
