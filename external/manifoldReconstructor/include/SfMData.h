/*
 * [Imported from external project]
 *
 * OpenMvgParser.cpp
 *
 *		Source: Manifold Reconstructor [https://github.com/andresax/Manifold-Reconstructor]
 *      Author: Andrea Romanoni
 */

#ifndef CAM_PARSERS_SFMPARSER_H_
#define CAM_PARSERS_SFMPARSER_H_

#include <glm.hpp>
#include "types_reconstructor.hpp"

struct SfMData {

  int numPoints_;
  int numCameras_;

  std::vector<glm::vec3> points_;
  std::vector<CameraType> camerasList_;
  std::vector<std::string> camerasPaths_;

  std::vector<std::vector<int> > camViewingPointN_;
  std::vector<std::vector<int> > pointsVisibleFromCamN_;
  std::vector<std::vector<glm::vec2> > point2DoncamViewingPoint_;

  int imageWidth_, imageHeight_;
};

#endif /* CAM_PARSERS_SFMPARSER_H_ */
