/*
 * [Imported from external project]
 *
 * OpenMvgParser.cpp
 *
 *		Source: Manifold Reconstructor [https://github.com/andresax/Manifold-Reconstructor]
 *      Author: Andrea Romanoni
 */

#ifndef CAM_PARSERS_OPENMVGPARSER_H_
#define CAM_PARSERS_OPENMVGPARSER_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <rapidjson/document.h>
#include "types_reconstructor.hpp"
#include "SfMData.h"


class OpenMvgParser {
public:
  OpenMvgParser(std::string path);
  virtual ~OpenMvgParser();
  virtual void parse();

  const SfMData& getSfmData() const {
    return sfm_data_;
  }

private:
  void parseViews(const std::map<int,glm::mat3> & intrinsics, const std::map<int,CameraType> & extrinsics, const std::map<int,int> &map_pos, const std::vector<int> &cam_ids);
  void parseIntrinsics(std::map<int,glm::mat3> & intrinsics);
  void parseExtrinsics(std::map<int,CameraType> & extrinsics, std::map<int,int> &map_pos, std::vector<int> &cam_ids);
  void parsePoints(const std::map<int,int> &map_pos);


  rapidjson::Document document_;
  std::string fileName_;
  std::ifstream fileStream_;
  SfMData sfm_data_;

  FILE* pFile_;
};

#endif /* CAM_PARSERS_OPENMVGPARSER_H_ */
