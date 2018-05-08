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

#include <getopt.h>
#include <stddef.h>
#include <iostream>

#include "SfMData.h"
#include "edge_graph_3d_utilities.hpp"
#include "output_point_cloud.hpp"

using namespace std;

/*
 * Converts a JSON in a point cloud in ply format. If original images are provided,
 * the point cloud will be colored.
 *
 * Input format:
 * ./json_to_ply [-i original_images/] input.json output.ply
 */
int main(int argc, char *argv[]) {
	char *input_json;
	char *output_ply;
	char *input_imgs = NULL;

	cout << "\n\njson_to_ply\n===========\n\nConverts a JSON in a point cloud in ply format.\n\n";

	bool invalid = false;

	  static const char *opt_string = "i:";
	  int opt = getopt(argc, argv, opt_string);
	  while (opt != -1) {
		switch (opt) {

		case 'i': input_imgs = optarg; break;

		default:
			invalid = true;
		}
		opt = getopt(argc, argv, opt_string);
	  }

	   if (argc - optind != 2) {
		invalid = true;
	  } else {
		  input_json = argv[optind];
		  output_ply = argv[optind+1];
	  }

	if(invalid) {
		cout << "Invalid input.\n\nCorrect format is:\n\n" << argv[0] << " [-i images/] input.json output.ply\n\tinput.json : input structure from motion data in JSON format\n\toutput.ply : output point cloud in PLY format\n\nOptional parameters:\n\n\t-i images/ : images originally used to produce the provided structure from motion data. If this argument is provided, json_to_ply will produce a colored point cloud\n\n";
		return 1;
	}

	cout << "Input: " << input_json << "\n";
	cout << "Output: " << output_ply << "\n";
	if(input_imgs != NULL)
		cout << "Original images: " << input_imgs << "\n";

	SfMData sfmd = read_sfm_data(input_json);
	if(input_imgs == NULL)
		output_point_cloud(sfmd, output_ply);
	else
		output_colored_point_cloud(sfmd, input_imgs, output_ply);

	return 0;
}


