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

#include <opencv2/core/types.hpp>
#include <cmath>
#include <limits>
#include <vector>

#include "polyline_graph_2d.hpp"
#include "polyline_graph_2d_hmap_impl.hpp"
#include "geometric_utilities.hpp"
#include "edge_graph_3d_utilities.hpp"
#include "filtering_close_plgps.hpp"
#include "glm.hpp"


float min2ddistsq_3dpoints(const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &a, const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &b) {
	float d=std::numeric_limits<float>::max();
	float cur_dist;

	for(int i=0; i < get<2>(a).size(); i++)
		for(int j=0; j < get<2>(b).size(); j++)
			if(get<2>(a)[i] == get<2>(b)[j] && get<1>(a)[i].polyline_id == get<1>(b)[j].polyline_id) {
				cur_dist = squared_2d_distance(get<1>(a)[i].plp.coords,get<1>(b)[j].plp.coords);
				if(cur_dist < d)
					d = cur_dist;
				break;
			}

	return d;
}

float closest2ddist_3dpoints(const vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &p3ds, const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &new_p3d) {
	float d=std::numeric_limits<float>::max();
	float cur_dist;
	ulong min_id=std::numeric_limits<ulong>::max();;

	for(ulong i=0; i < p3ds.size(); i++) {
		cur_dist = min2ddistsq_3dpoints(p3ds[i],new_p3d);
		if(cur_dist < d) {
			d = cur_dist;
			min_id = i;
		}
	}

	return sqrt(d);
}

#define CELLSIZE 3

bool getcellvalue(bool ** bmap, const glm::vec2 &v)
{
	return bmap[int(v.y/CELLSIZE)][int(v.x/CELLSIZE)];
}

void setcellvalue(bool ** bmap, const glm::vec2 &v, const bool newval)
{
	bmap[int(v.y/CELLSIZE)][int(v.x/CELLSIZE)] = newval;
}

bool is_new_point(const vector<bool**> &bmaps, const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &new_3dpt) {
	for(int i=0; i < get<1>(new_3dpt).size(); i++)
		if(!getcellvalue(bmaps[get<2>(new_3dpt)[i]],get<1>(new_3dpt)[i].plp.coords))
			return true;
	return false;
}

void set_new_point(vector<bool**> &bmaps, const std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>> &new_3dpt) {
	for(int i=0; i < get<1>(new_3dpt).size(); i++)
		setcellvalue(bmaps[get<2>(new_3dpt)[i]],get<1>(new_3dpt)[i].plp.coords,true);
}

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> filter_3d_points_close_2d_array(const vector<PolyLineGraph2DHMapImpl> &plgs, const Size &imgsz, const vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &p3ds) {
	Size sz(int(ceil((float)imgsz.width/CELLSIZE)),int(ceil((float)imgsz.height/CELLSIZE)));

	vector<bool**> bmaps(plgs.size());

	for(int i=0; i < plgs.size(); i++) {
		bmaps[i] = create_2D_array<bool>(sz.height,sz.width);
		for(int row=0; row < sz.height; row++)
			for(int col=0; col < sz.width; col++)
				bmaps[i][row][col] = false;
	}

	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> res;

	for(const auto &p3d : p3ds)
		if(is_new_point(bmaps,p3d)){
			set_new_point(bmaps,p3d);
			res.push_back(p3d);
		}


	for(int i=0; i < plgs.size(); i++)
		delete_2D_array(bmaps[i]);

	return res;
}

vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> filter3dpoints_close2d(const vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> &p3ds) {
	vector<std::tuple<glm::vec3, vector<PolyLineGraph2D::plg_point>, vector<int>>> res;

	for(const auto &p3d : p3ds)
		if(closest2ddist_3dpoints(res,p3d) >= 3.0)
			res.push_back(p3d);

	return res;
}



