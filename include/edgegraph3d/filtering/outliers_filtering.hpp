/*
 * outliers_filtering.hpp
 *
 *  Created on: Apr 14, 2018
 *      Author: root
 */

#ifndef INCLUDE_EDGEGRAPH3D_FILTERING_OUTLIERS_FILTERING_HPP_
#define INCLUDE_EDGEGRAPH3D_FILTERING_OUTLIERS_FILTERING_HPP_

#include "global_defines.hpp"
#include "global_switches.hpp"

struct SfMData;

#define FILTER_3VIEWS_AMOUNT 3

void filter(SfMData &sfmd, const int first_edgepoint);
void filter(SfMData &sfmd, const int first_edgepoint, const float gn_max_mse);
void filter(SfMData &sfmd, const int first_edgepoint, const int forced_min_filter);
void filter(SfMData &sfmd, const int first_edgepoint, const float gn_max_mse, const int forced_min_filter);

#endif /* INCLUDE_EDGEGRAPH3D_FILTERING_OUTLIERS_FILTERING_HPP_ */
