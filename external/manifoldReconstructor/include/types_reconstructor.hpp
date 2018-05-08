//  Copyright 2014 Andrea Romanoni
//
//  This file is part of edgePointSpaceCarver.
//
//  edgePointSpaceCarver is free software: you can redistribute it
//  and/or modify it under the terms of the GNU General Public License as
//  published by the Free Software Foundation, either version 3 of the
//  License, or (at your option) any later version see
//  <http://www.gnu.org/licenses/>..
//
//  edgePointSpaceCarver is distributed in the hope that it will be
//  useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

/*
 * [Imported from external project]
 *
 * OpenMvgParser.cpp
 *
 *		Source: Manifold Reconstructor [https://github.com/andresax/Manifold-Reconstructor]
 *      Author: Andrea Romanoni
 */

/**
* Header-only file with various types, especially those related to the CGAL library
*/
#ifndef TYPES_RECONSTR_HPP_
#define TYPES_RECONSTR_HPP_

#include <string>
#include <glm.hpp>
#include <Eigen/Core>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/intersections.h>

#include <CGAL/algorithm.h>

#include "Delaunay3DCellInfo.h"
#include "Delaunay3DVertexInfo.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;
typedef K::Segment_3 Segment;

typedef CGAL::Triangulation_vertex_base_with_info_3<Delaunay3DVertexInfo, K> Vb;

typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vb> Vbh;
typedef CGAL::Triangulation_cell_base_with_info_3<Delaunay3DCellInfo, K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vbh, Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Dt;
typedef CGAL::Triangulation_hierarchy_3<Dt> Delaunay3;
typedef Delaunay3::Point PointD3;
typedef Delaunay3::Vertex_handle Vertex3D_handle;

typedef std::pair<PointD3, float> DistanceWeight;

typedef std::set<Delaunay3DCellInfo::FSConstraint, Delaunay3DCellInfo::LtFSConstraint> SetConstraints;

struct CameraType {
  glm::mat3 intrinsics;
  glm::mat3 rotation;
  glm::vec3 translation;
  glm::mat4 cameraMatrix;
  glm::vec3 center;
  glm::mat4 mvp;

  std::string pathImage;

  int imageWidth;
  int imageHeight;

  std::vector<int> visiblePoints;
};

struct PointParser {
  float x;
  float y;
  float z;

  int R;
  int G;
  int B;

  //position of the feature in the corresponding image;
  //the center of the image plane is the origin
  std::vector<float> viewingCamerasX;
  std::vector<float> viewingCamerasY;

  std::vector<int> viewingCamerasIndices;
};

struct sortTetByIntersection {
  inline bool operator()(const Delaunay3::Cell_handle& i, const Delaunay3::Cell_handle& j) {
   return i->info().getVoteCountProb() < j->info().getVoteCountProb();
  }
};


struct PointReconstruction {
  PointD3 position;
  Vertex3D_handle vertexHandle;
  int idVertex;
  bool new_;
  std::vector<int> viewingCams;

  PointReconstruction() {
    position = PointD3(0.0, 0.0, 0.0);
    new_ = true;
    idVertex = -1;
  }
};

struct CamReconstruction {
  PointD3 position;
  Vertex3D_handle vertexHandle;
  std::vector<int> visiblePoints;
  std::vector<int> newVisiblePoints;

  CamReconstruction() {
    position = PointD3(0.0, 0.0, 0.0);
  }
};



#endif /* TYPES_HPP_ */
