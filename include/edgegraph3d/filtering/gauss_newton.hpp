/*
 * gauss_newton.hpp
 *
 *  Created on: Apr 14, 2018
 *      Author: root
 */

#ifndef INCLUDE_EDGEGRAPH3D_FILTERING_GAUSS_NEWTON_HPP_
#define INCLUDE_EDGEGRAPH3D_FILTERING_GAUSS_NEWTON_HPP_

#include <vector>

#include "global_defines.hpp"
#include "global_switches.hpp"

struct SfMData;

#define GN_MAX_MSE 2.25

void gaussNewtonFiltering(SfMData &sfm_data_, std::vector<bool>& inliers, const float gn_max_mse);

#endif /* INCLUDE_EDGEGRAPH3D_FILTERING_GAUSS_NEWTON_HPP_ */
