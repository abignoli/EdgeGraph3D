
#include "outliers_filtering.hpp"

#include <iostream>
#include <vector>

#include "SfMData.h"
#include "gauss_newton.hpp"

using namespace std;

#define INVALID_FORCED_MIN_FILTER -1

void compute_ray_stats(SfMData &sfmd, const std::vector<bool> &inliers, float &average_rays, int &median_ray_amount) {
	std::vector<int> point_rays_amount_distribution(sfmd.camerasList_.size(),0);

	int count = 0;
	unsigned long amount_of_rays=0;
	for(int i=0; i < sfmd.numPoints_; i++)
		if(inliers[i]) {
			count++;
			point_rays_amount_distribution[sfmd.camViewingPointN_[i].size()-1]++;
			amount_of_rays += sfmd.camViewingPointN_[i].size();
		}

	average_rays = amount_of_rays / ((float) count);

	int m_amount_of_rays=0;
	for(median_ray_amount=0; median_ray_amount < sfmd.numCameras_; median_ray_amount++) {
		m_amount_of_rays += point_rays_amount_distribution[median_ray_amount];
		if(m_amount_of_rays >= count/2)
			break;
	}

}

vector<bool> compute_inliers(SfMData &sfm_data_, const int first_edgepoint, const float gn_max_mse, const int forced_min_filter) {
	vector<bool> inliers;
	gaussNewtonFiltering(sfm_data_, inliers, gn_max_mse);

	float average_rays;
	int median_ray_amount;
	compute_ray_stats(sfm_data_, inliers, average_rays,median_ray_amount);

	int gn_count = 0;
	for (int curpt = 0; curpt < sfm_data_.points_.size(); curpt++)
		if (!inliers[curpt])
			gn_count++;

	std::cout << "Gauss-Newton filtered: " << gn_count << std::endl;

	int intended_view_filter= (FILTER_3VIEWS_AMOUNT >= median_ray_amount/2 - 1) ? FILTER_3VIEWS_AMOUNT : (median_ray_amount/2 - 1);

	if(forced_min_filter > INVALID_FORCED_MIN_FILTER)
		intended_view_filter = forced_min_filter;

    std::cout << "Filtering points with ID >= " << first_edgepoint << " with less than " << intended_view_filter+1 << " observations\n";

	for (int curpt = first_edgepoint;
			curpt < sfm_data_.points_.size(); curpt++)
		inliers[curpt] = inliers[curpt] && (sfm_data_.camViewingPointN_[curpt].size() > intended_view_filter);

	return inliers;
}

void removeOutliers(SfMData &sfmd, const vector<bool> &inliers) {
	SfMData res;

	res.camerasList_ = sfmd.camerasList_;
	res.camerasPaths_ = sfmd.camerasPaths_;
	res.numCameras_ = sfmd.numCameras_;

	res.imageWidth_ = sfmd.imageWidth_;
	res.imageHeight_ = sfmd.imageHeight_;

	res.pointsVisibleFromCamN_.resize(res.numCameras_);

	int curpt = 0;
	for(int i=0; i < inliers.size(); i++)
		if(inliers[i]) {
			res.points_.push_back(sfmd.points_[i]);
			res.camViewingPointN_.push_back(sfmd.camViewingPointN_[i]);
			res.point2DoncamViewingPoint_.push_back(sfmd.point2DoncamViewingPoint_[i]);
			for(const auto cam_id : sfmd.camViewingPointN_[i])
				res.pointsVisibleFromCamN_[cam_id].push_back(curpt);
			curpt++;
		}

	res.numPoints_ = res.points_.size();

	sfmd = res;
}

void filter(SfMData &sfmd, const int first_edgepoint) {
	filter(sfmd, first_edgepoint, GN_MAX_MSE, INVALID_FORCED_MIN_FILTER);
}

void filter(SfMData &sfmd, const int first_edgepoint, const float gn_max_mse) {
	return filter(sfmd, first_edgepoint, gn_max_mse, INVALID_FORCED_MIN_FILTER);
}

void filter(SfMData &sfmd, const int first_edgepoint, const int forced_min_filter) {
	return filter(sfmd, first_edgepoint, GN_MAX_MSE, forced_min_filter);
}

void filter(SfMData &sfmd, const int first_edgepoint, const float gn_max_mse, const int forced_min_filter) {
	cout << "Filtering... ";
	int init_points = sfmd.points_.size();
	auto inliers = compute_inliers(sfmd, first_edgepoint, gn_max_mse, forced_min_filter);
	removeOutliers(sfmd, inliers);
	int final_points = sfmd.points_.size();
	cout << "Removed " << init_points-final_points << " points.\n\n";
	cout << "Final amount of computed 3D points: " << sfmd.points_.size() << "\n";
}
