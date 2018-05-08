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


#include "polyline_graph_2d_hmap_impl.hpp"

#include <cmath>
#include <set>
#include <utility>

#include "polyline_graph_2d_hmap_impl.hpp"
#include "polyline_2d_map.hpp"

PolyLine2DMap::PolyLine2DMap(const PolyLineGraph2DHMapImpl &plg, const Size img_sz, const float cell_dim) : plg(plg), img_sz(img_sz), cell_dim(cell_dim) {
	mapsz = Size(int(ceil(img_sz.width/cell_dim)),int(ceil(img_sz.height/cell_dim)));
	pls_id_maps = create_2D_array<vector<ulong>>(mapsz.height,mapsz.width);
	for(ulong row=0; row < mapsz.height; row++)
		for(ulong col=0; col < mapsz.width; col++)
			pls_id_maps[row][col]=vector<ulong>();
	for(ulong pl_id=0; pl_id < plg.polylines.size(); pl_id++)
		add_polyline(pl_id);
}

PolyLine2DMap::~PolyLine2DMap() {}

void PolyLine2DMap::add_polyline(const ulong pl_id) {
	if(!(plg.is_valid_polyline(pl_id)))
		return;
	set<pair<ulong,ulong>> intersectedcells = plg.polylines[pl_id].get_intersectedcells_2dmap_set(cell_dim);
	for(const auto &intersectedcell : intersectedcells)
		pls_id_maps[intersectedcell.second][intersectedcell.first].push_back(pl_id);
}
