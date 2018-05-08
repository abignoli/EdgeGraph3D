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

#include <iostream>
#include <string>

#include "SfMData.h"
#include "edge_graph_3d_utilities.hpp"
#include "output_sfm_data.hpp"
#include "transform_coordinate_system.hpp"

using namespace std;

/*
 * Given a SfMData JSON and the list of true camera poses (or the camera poses in the target
 * coordinate system), returns a JSON adapted to the true (target) coordinate system. Useful
 * for comparisons with ground truth. Camera poses should be expressed in a ASCII file
 * containing on each line the triplet x y z of coordinates of a camera. The order of the
 * cameras need to match the one in the input JSON.
 *
 * Input format:
 * ./coordinate_system_transform input.json camera_poses.txt output.json
 */
int main(int argc, char *argv[]) {
	cout << "\n\ncoordinate_system_transform\n===========================\n\nTransform the coordinate system in the input JSON to match the provided camera poses.\n\n";

	if(argc != 4) {
		cout << "Invalid input.\n\nCorrect format is:\n\n" << argv[0] << " input.json camera_poses.txt output.json\n\tinput.json : input structure from motion data in JSON format\n\tcamera_poses.txt : files containing the target camera poses, one per line, in the order defined by input.json, in the format \"x y z\"\n\toutput.json : output transformed coordinate system structure from motion data in JSON format\n\n";
		return 1;
	}

	char *input_json = argv[1];
	char *target_camera_poses = argv[2];
	char *output_json = argv[3];

	cout << "Input: " << input_json << "\n";
	cout << "Target camera poses: " << target_camera_poses << "\n";
	cout << "Output: " << output_json << "\n";

	SfMData sfmd = read_sfm_data(input_json);
	update_sfmd_from_real_camera_positions(sfmd, string(target_camera_poses));
	output_sfm_data(input_json, sfmd, output_json);
}


