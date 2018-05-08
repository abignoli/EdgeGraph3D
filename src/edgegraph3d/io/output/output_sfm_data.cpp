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


#include "output_sfm_data.hpp"

#include <sys/types.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../../../../external/manifoldReconstructor/include/Exceptions.hpp"
#include "SfMData.h"
#include "types_reconstructor.hpp"
#include "glm.hpp"

#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"
#include "rapidjson/ostreamwrapper.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/rapidjson.h"


using namespace std;
using namespace rapidjson;

rapidjson::Value get_sfmd_version_value(const rapidjson::Document &sfmd_doc, rapidjson::Document::AllocatorType &a) {
	if (!sfmd_doc.HasMember("sfm_data_version"))
	      throw JsonAccessException("JsonAccessException--> error while querying HasMember(sfm_data_version)");
	return rapidjson::Value(sfmd_doc["sfm_data_version"],a) ;
}

rapidjson::Value get_sfmd_root_path_value(const rapidjson::Document &sfmd_doc, rapidjson::Document::AllocatorType &a) {
	if (!sfmd_doc.HasMember("root_path"))
	      throw JsonAccessException("JsonAccessException--> error while querying HasMember(root_path)");
	return rapidjson::Value(sfmd_doc["root_path"],a);
}

rapidjson::Value get_sfmd_intrinsics_value(const rapidjson::Document &sfmd_doc, rapidjson::Document::AllocatorType &a) {
	if (!sfmd_doc.HasMember("intrinsics"))
	      throw JsonAccessException("JsonAccessException--> error while querying HasMember(intrinsics)");
	return rapidjson::Value(sfmd_doc["intrinsics"],a);
}

rapidjson::Value get_sfmd_views_value(const rapidjson::Document &sfmd_doc, rapidjson::Document::AllocatorType &a) {
	if (!sfmd_doc.HasMember("views"))
	      throw JsonAccessException("JsonAccessException--> error while querying HasMember(views)");
	return rapidjson::Value(sfmd_doc["views"],a);
}

rapidjson::Value get_sfmd_control_points_value(const rapidjson::Document &sfmd_doc, rapidjson::Document::AllocatorType &a) {
	if (!sfmd_doc.HasMember("control_points"))
	      throw JsonAccessException("JsonAccessException--> error while querying HasMember(control_points)");
	return rapidjson::Value(sfmd_doc["control_points"],a);
}

rapidjson::Value convert_glmmat3_to_rapidjson_val(const glm::mat3 m, rapidjson::Document::AllocatorType &a) {
	rapidjson::Value rjs_m(rapidjson::kArrayType);
	for(int row=0; row<3; row++) {
		rapidjson::Value row_vals(rapidjson::kArrayType);
		for(int col=0; col<3;col++)
			row_vals.PushBack(rapidjson::Value(m[row][col]),a);
		rjs_m.PushBack(row_vals,a);
	}
	return rjs_m;
}

rapidjson::Value convert_glmvec3_to_rapidjson_val(const glm::vec3 &v, rapidjson::Document::AllocatorType &a) {
	rapidjson::Value center(rapidjson::kArrayType);
	for(int row=0; row<3; row++)
		center.PushBack(rapidjson::Value(v[row]),a);
	return center;
}

rapidjson::Value convert_glmvec2_to_rapidjson_val(const glm::vec2 &v, rapidjson::Document::AllocatorType &a) {
	rapidjson::Value center(rapidjson::kArrayType);
	center.PushBack(rapidjson::Value(v.x),a);
	center.PushBack(rapidjson::Value(v.y),a);
	return center;
}

rapidjson::Value get_cam_extrinsics(const CameraType &c, int key_id, rapidjson::Document::AllocatorType &a) {
	rapidjson::Value extrinsics;
	extrinsics.SetObject();

	extrinsics.AddMember("key",key_id,a);

	rapidjson::Value pose;
	pose.SetObject();

	rapidjson::Value rotation = convert_glmmat3_to_rapidjson_val(c.rotation,a);
	pose.AddMember("rotation",rotation,a);

	rapidjson::Value center = convert_glmvec3_to_rapidjson_val(c.center,a);
	pose.AddMember("center",center,a);

	extrinsics.AddMember("value",pose,a);

	return extrinsics;
}

rapidjson::Value get_extrinsics(const SfMData &sfmd, rapidjson::Document::AllocatorType &a) {
	rapidjson::Value extrinsics(rapidjson::kArrayType);
	for(int i=0; i < sfmd.numCameras_; i++)
		extrinsics.PushBack(get_cam_extrinsics(sfmd.camerasList_[i],i,a),a);
	return extrinsics;
}

#define OUTPUT_SFMD_FEATURE_ID 0
rapidjson::Value get_refpoint(const SfMData &sfmd, const ulong refpoint_id, rapidjson::Document::AllocatorType &a) {
	rapidjson::Value refpoint;
	refpoint.SetObject();

	refpoint.AddMember("key",refpoint_id,a);

	rapidjson::Value refpoint_value;
	refpoint_value.SetObject();

	refpoint_value.AddMember("X",convert_glmvec3_to_rapidjson_val(sfmd.points_[refpoint_id],a),a);

	rapidjson::Value observations(rapidjson::kArrayType);
	for(int i=0; i < sfmd.camViewingPointN_[refpoint_id].size();i++) {
		rapidjson::Value observation;
		observation.SetObject();
		const int cam_id = sfmd.camViewingPointN_[refpoint_id][i];
		const glm::vec2 &coords = sfmd.point2DoncamViewingPoint_[refpoint_id][i];
		observation.AddMember("key",cam_id,a);

		rapidjson::Value observation_value;
		observation_value.SetObject();
		observation_value.AddMember("id_feat",OUTPUT_SFMD_FEATURE_ID,a);
		observation_value.AddMember("x",convert_glmvec2_to_rapidjson_val(coords,a),a);
		observation.AddMember("value",observation_value,a);

		observations.PushBack(observation,a);
	}
	refpoint_value.AddMember("observations",observations,a);

	refpoint.AddMember("value",refpoint_value,a);

	return refpoint;
}

rapidjson::Value get_structure(const SfMData &sfmd, rapidjson::Document::AllocatorType &a) {
	rapidjson::Value structure(rapidjson::kArrayType);
	for(ulong refpoint_id=0; refpoint_id < sfmd.numPoints_; refpoint_id++)
		structure.PushBack(get_refpoint(sfmd,refpoint_id,a),a);
	return structure;
}

rapidjson::Value get_structure(const SfMData &sfmd, rapidjson::Document::AllocatorType &a, const vector<bool> &inliers) {
	rapidjson::Value structure(rapidjson::kArrayType);
	for(ulong refpoint_id=0; refpoint_id < sfmd.numPoints_; refpoint_id++)
		if(inliers[refpoint_id])
			structure.PushBack(get_refpoint(sfmd,refpoint_id,a),a);
	return structure;
}

void output_sfm_data(const char *original_sfmd_filename, const SfMData &sfmd, const std::string &output_file) {
	ifstream original_sfmd_file;
	original_sfmd_file.open(original_sfmd_filename);

	rapidjson::Document original_sfmd_doc;

	rapidjson::IStreamWrapper isw(original_sfmd_file);
	original_sfmd_doc.ParseStream(isw);

	  try {
	    if (!original_sfmd_doc.IsObject())
	      throw JsonParseException("JsonParseException--> the json file " + std::string(original_sfmd_filename) + " is not valid");
	  } catch (JsonParseException& e) {
	    std::cerr << e.what() << std::endl;
	    std::cout << e.what() << std::endl;
	  }

	rapidjson::Document jsonDoc;
	jsonDoc.SetObject();

	rapidjson::Document::AllocatorType& allocator = jsonDoc.GetAllocator();


	rapidjson::Value sfmdValue;
	sfmdValue.SetObject();

	sfmdValue.AddMember("sfm_data_version",get_sfmd_version_value(original_sfmd_doc,allocator),allocator);
	sfmdValue.AddMember("root_path",get_sfmd_root_path_value(original_sfmd_doc,allocator),allocator);
	sfmdValue.AddMember("views",get_sfmd_views_value(original_sfmd_doc,allocator),allocator);
	sfmdValue.AddMember("intrinsics",get_sfmd_intrinsics_value(original_sfmd_doc,allocator),allocator);
	sfmdValue.AddMember("extrinsics",get_extrinsics(sfmd,allocator),allocator);
	sfmdValue.AddMember("structure",get_structure(sfmd,allocator),allocator);

	// control points
	sfmdValue.AddMember("control_points",get_sfmd_control_points_value(original_sfmd_doc,allocator),allocator);



	ofstream ofs(output_file);
	OStreamWrapper osw(ofs);
	PrettyWriter<OStreamWrapper> writer(osw);
	sfmdValue.Accept(writer);
	ofs.close();
}



void output_filtered_sfm_data(const char *original_sfmd_filename, const SfMData &sfmd, const vector<bool> &inliers, const std::string &output_file) {
	ifstream original_sfmd_file;
	original_sfmd_file.open(original_sfmd_filename);

	rapidjson::Document original_sfmd_doc;

	rapidjson::IStreamWrapper isw(original_sfmd_file);
	original_sfmd_doc.ParseStream(isw);

	  try {
		if (!original_sfmd_doc.IsObject())
		  throw JsonParseException("JsonParseException--> the json file " + std::string(original_sfmd_filename) + " is not valid");
	  } catch (JsonParseException& e) {
		std::cerr << e.what() << std::endl;
		std::cout << e.what() << std::endl;
	  }

	rapidjson::Document jsonDoc;
	jsonDoc.SetObject();

	rapidjson::Document::AllocatorType& allocator = jsonDoc.GetAllocator();


	rapidjson::Value sfmdValue;
	sfmdValue.SetObject();

	sfmdValue.AddMember("sfm_data_version",get_sfmd_version_value(original_sfmd_doc,allocator),allocator);
	sfmdValue.AddMember("root_path",get_sfmd_root_path_value(original_sfmd_doc,allocator),allocator);
	sfmdValue.AddMember("views",get_sfmd_views_value(original_sfmd_doc,allocator),allocator);
	sfmdValue.AddMember("intrinsics",get_sfmd_intrinsics_value(original_sfmd_doc,allocator),allocator);
	sfmdValue.AddMember("extrinsics",get_extrinsics(sfmd,allocator),allocator);
	sfmdValue.AddMember("structure",get_structure(sfmd,allocator,inliers),allocator);

	// control points
	sfmdValue.AddMember("control_points",get_sfmd_control_points_value(original_sfmd_doc,allocator),allocator);

	ofstream ofs(output_file);
	OStreamWrapper osw(ofs);
	PrettyWriter<OStreamWrapper> writer(osw);
	sfmdValue.Accept(writer);
	ofs.close();
	original_sfmd_file.close();
}

