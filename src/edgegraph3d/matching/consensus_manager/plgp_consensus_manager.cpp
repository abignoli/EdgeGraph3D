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

#include "SfMData.h"
#include "plgp_consensus_manager.hpp"


using namespace cv;

PLGPConsensusManager::PLGPConsensusManager(vector<Mat> &input_imgs, const SfMData &input_sfmd, const Size &input_sz, const Mat** all_fundamental_matrices, vector<PolyLineGraph2DHMapImpl> &plgs) : imgs(input_imgs), sfmd(input_sfmd),img_size(input_sz), all_fundamental_matrices(all_fundamental_matrices), plgs(plgs) {}


