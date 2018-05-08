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


#ifndef INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_POLYLINE_2D_MAP_HPP_
#define INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_POLYLINE_2D_MAP_HPP_

#include <opencv2/core/types.hpp>
#include <vector>

#include "edge_graph_3d_utilities.hpp"
#include "global_defines.hpp"
#include "global_switches.hpp"

class PolyLineGraph2DHMapImpl;

using namespace std;
using namespace cv;

class PolyLine2DMap {
private:
	void add_polyline(const ulong pl_id);
public:
	const PolyLineGraph2DHMapImpl &plg;
	vector<ulong> **pls_id_maps;
	const float cell_dim;
	const Size img_sz;
	const Size mapsz;

	PolyLine2DMap(const PolyLineGraph2DHMapImpl &plg, const Size img_sz, const float cell_dim);
	~PolyLine2DMap();
};

#endif /* INCLUDE_EDGEGRAPH3D_MATCHING_PLG_MATCHING_POLYLINE_2D_MAP_HPP_ */
