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


#include "test_utilities.hpp"

#include <iostream>
#include <vector>

#include "SfMData.h"
#include "geometric_utilities.hpp"
#include "glm.hpp"

void test_reprojection(const SfMData &sfmd, const int cam_id, const glm::vec3 &p3d, const glm::vec2 &prj) {
	glm::vec2 computed_prj = compute_projection(sfmd,cam_id,p3d);
	cout << "Camera " << cam_id << ": saved_projection = (" << prj[0] << "," << prj[1] << "), computed projection (" << computed_prj[0] << "," << computed_prj[1] << ")" << endl;
}

void test_refpoint_projections(const SfMData &sfmd,const int refpoint) {
	cout << "Testing reprojections of refpoint " << refpoint << endl;
	const glm::vec3 &p3d = sfmd.points_[refpoint];
	for(int i=0; i < sfmd.camViewingPointN_[refpoint].size(); i++) {
		const int cam_id = sfmd.camViewingPointN_[refpoint][i];
		const glm::vec2 &coords = sfmd.point2DoncamViewingPoint_[refpoint][i];
		test_reprojection(sfmd,cam_id,p3d,coords);
	}
}

