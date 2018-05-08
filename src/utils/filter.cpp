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
#include <cstdlib>
#include <iostream>

#include "SfMData.h"
#include "edge_graph_3d_utilities.hpp"
#include "outliers_filtering.hpp"
#include "output_sfm_data.hpp"

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
	int starting_point = 0;
	float gn_max_mse = 0;
	int force_views_filtering = -1;
	char *output_json;

	bool invalid = false;

	cout << "\n\nfilter\n======\n\nFilters outliers in structure from motion data.\n\n";

	  static const char *opt_string = "s:e:f:";
	  int opt = getopt(argc, argv, opt_string);
	  while (opt != -1) {
		switch (opt) {

		case 's': starting_point = atoi(optarg);
			if(starting_point < 0)
				cout << "starting_point_id must be positive or null.\n";
			break;

		case 'e': gn_max_mse = atof(optarg);
			if(gn_max_mse < 0)
				cout << "Gauss Netwon maximum mean squared error must be positive.\n";
			break;

		case 'f': force_views_filtering = atoi(optarg);
			if(force_views_filtering < 0)
				cout << "force_views_filtering must be positive or null.\n";
			break;


		default:
			invalid = true;
		}
		opt = getopt(argc, argv, opt_string);
	  }

	   if (argc - optind != 2) {
		invalid = true;
	  } else {
		  input_json = argv[optind];
		  output_json = argv[optind+1];
	  }

	if(invalid) {
		cout << "Invalid input.\n\nCorrect format is:\n\n" << argv[0] << "[-s starting_point_id | -e gn_max_mse | -f force_views_filtering] input.json output.json\n\tinput.json : input structure from motion data in JSON format\n\toutput.json : output inliers point cloud in JSON format\n\nOptional parameters:\n\n\t-s starting_point_id : filters points starting from starting_point_id (default is zero) based on number of observations\n\t-e gn_max_mse : maximum mean squared error to be tollerated by Gauss Newton procedure in outliers filtering\n\t-f force_views_filtering : removes all points after starting_point_id with a count of observations smaller than or equal to force_views_filtering\n\n";
		return 1;
	}

	cout << "Input: " << input_json << "\n";
	cout << "Output: " << output_json << "\n\n";

	SfMData sfmd = read_sfm_data(input_json);
	if(force_views_filtering < 0) {
		if(gn_max_mse <= 0)
			filter(sfmd, starting_point);
		else
			filter(sfmd, starting_point, gn_max_mse);
	} else {
		if(gn_max_mse <= 0)
			filter(sfmd, starting_point, force_views_filtering);
		else
			filter(sfmd, starting_point, gn_max_mse, force_views_filtering);
	}
	output_sfm_data(input_json, sfmd, output_json);
	return 0;
}
