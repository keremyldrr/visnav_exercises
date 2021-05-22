/**
BSD 3-Clause License

Copyright (c) 2018, Vladyslav Usenko and Nikolaus Demmel.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

#include <fstream>

#include <tbb/task_scheduler_init.h>

#include <ceres/ceres.h>

#include <opengv/absolute_pose/CentralAbsoluteAdapter.hpp>
#include <opengv/absolute_pose/methods.hpp>
#include <opengv/relative_pose/CentralRelativeAdapter.hpp>
#include <opengv/sac/Ransac.hpp>
#include <opengv/sac_problems/absolute_pose/AbsolutePoseSacProblem.hpp>
#include <opengv/triangulation/methods.hpp>

#include <visnav/common_types.h>
#include <visnav/serialization.h>

#include <visnav/reprojection.h>
#include <visnav/local_parameterization_se3.hpp>

#include <visnav/tracks.h>
using namespace opengv;

namespace visnav {

// save map with all features and matches
void save_map_file(const std::string& map_path, const Corners& feature_corners,
                   const Matches& feature_matches,
                   const FeatureTracks& feature_tracks,
                   const FeatureTracks& outlier_tracks, const Cameras& cameras,
                   const Landmarks& landmarks) {
  {
    std::ofstream os(map_path, std::ios::binary);

    if (os.is_open()) {
      cereal::BinaryOutputArchive archive(os);
      archive(feature_corners);
      archive(feature_matches);
      archive(feature_tracks);
      archive(outlier_tracks);
      archive(cameras);
      archive(landmarks);

      size_t num_obs = 0;
      for (const auto& kv : landmarks) {
        num_obs += kv.second.obs.size();
      }
      std::cout << "Saved map as " << map_path << " (" << cameras.size()
                << " cameras, " << landmarks.size() << " landmarks, " << num_obs
                << " observations)" << std::endl;
    } else {
      std::cout << "Failed to save map as " << map_path << std::endl;
    }
  }
}

// load map with all features and matches
void load_map_file(const std::string& map_path, Corners& feature_corners,
                   Matches& feature_matches, FeatureTracks& feature_tracks,
                   FeatureTracks& outlier_tracks, Cameras& cameras,
                   Landmarks& landmarks) {
  {
    std::ifstream is(map_path, std::ios::binary);

    if (is.is_open()) {
      cereal::BinaryInputArchive archive(is);
      archive(feature_corners);
      archive(feature_matches);
      archive(feature_tracks);
      archive(outlier_tracks);
      archive(cameras);
      archive(landmarks);

      size_t num_obs = 0;
      for (const auto& kv : landmarks) {
        num_obs += kv.second.obs.size();
      }
      std::cout << "Loaded map from " << map_path << " (" << cameras.size()
                << " cameras, " << landmarks.size() << " landmarks, " << num_obs
                << " observations)" << std::endl;
    } else {
      std::cout << "Failed to load map from " << map_path << std::endl;
    }
  }
}

// Create new landmarks from shared feature tracks if they don't already exist.
// The two cameras must be in the map already.
// Returns the number of newly created landmarks.
int add_new_landmarks_between_cams(const FrameCamId& fcid0,
                                   const FrameCamId& fcid1,
                                   const Calibration& calib_cam,
                                   const Corners& feature_corners,
                                   const FeatureTracks& feature_tracks,
                                   const Cameras& cameras,
                                   Landmarks& landmarks) {
  // shared_track_ids will contain all track ids shared between the two images,
  // including existing landmarks
  std::vector<TrackId> shared_track_ids;

  // find shared feature tracks
  const std::set<FrameCamId> fcids = {fcid0, fcid1};
  if (!GetTracksInImages(fcids, feature_tracks, shared_track_ids)) {
    return 0;
  }

  bearingVectors_t vec1;

  bearingVectors_t vec2;
  // run method 1
  std::vector<TrackId> new_track_ids;

  for (unsigned long i = 0; i < shared_track_ids.size(); i++) {
    if (landmarks.find(shared_track_ids[i]) != landmarks.end()) continue;

    FeatureTrack track = feature_tracks.at(shared_track_ids[i]);
    auto f0 = track[fcid0];
    auto f1 = track[fcid1];
    auto c0 = feature_corners.at(fcid0).corners[f0];
    auto c1 = feature_corners.at(fcid1).corners[f1];
    bearingVector_t p0(calib_cam.intrinsics[fcid0.cam_id]->unproject(c0));
    bearingVector_t p1(calib_cam.intrinsics[fcid1.cam_id]->unproject(c1));
    vec1.push_back(p0);
    vec2.push_back(p1);
    new_track_ids.push_back(shared_track_ids[i]);

    /// triangulate and add to landmarks here
  }

  auto fc1to0 = cameras.at(fcid0).T_w_c.inverse() * cameras.at(fcid1).T_w_c;
  relative_pose::CentralRelativeAdapter adapter(
      vec1, vec2, fc1to0.translation(), fc1to0.rotationMatrix());

  for (unsigned long j = 0; j < new_track_ids.size(); j++) {
    point_t point = triangulation::triangulate(adapter, j);
    Landmark temp;
    temp.p = cameras.at(fcid0).T_w_c * point;
    FeatureTrack track = feature_tracks.at(new_track_ids[j]);
    // temp.obs.insert(std::make_pair(fcid0, track.at(fcid0)));

    for (auto kv : cameras) {
      auto cid = kv.first;

      if (track.find(cid) != track.end())  // if not in
      {
        // temp.obs.insert(std::make_pair(kv.first, track.at(cid)));
        temp.obs[cid] = track.at(cid);
      }
    }
    landmarks[new_track_ids[j]] = temp;
  }
  // TODO SHEET 4: Triangulate all new features and add to the map

  return new_track_ids.size();
}

// Initialize the scene from a stereo pair, using the known transformation from
// camera calibration. This adds the inital two cameras and triangulates shared
// landmarks.
// Note: in principle we could also initialize a map from another images pair
// using the transformation from the pairwise matching with the 5-point
// algorithm. However, using a stereo pair has the advantage that the map is
// initialized with metric scale.
bool initialize_scene_from_stereo_pair(const FrameCamId& fcid0,
                                       const FrameCamId& fcid1,
                                       const Calibration& calib_cam,
                                       const Corners& feature_corners,
                                       const FeatureTracks& feature_tracks,
                                       Cameras& cameras, Landmarks& landmarks) {
  // check that the two image ids refer to a stereo pair
  if (!(fcid0.frame_id == fcid1.frame_id && fcid0.cam_id != fcid1.cam_id)) {
    std::cerr << "Images " << fcid0 << " and " << fcid1
              << " don't form a stereo pair. Cannot initialize." << std::endl;
    return false;
  }

  // TODO SHEET 4: Initialize scene (add initial cameras and landmarks)

  cameras[fcid0].T_w_c = Sophus::SE3d(Eigen::Matrix4d::Identity());
  cameras[fcid1].T_w_c = calib_cam.T_i_c[fcid1.cam_id];  // TODO : Make sure

  add_new_landmarks_between_cams(fcid0, fcid1, calib_cam, feature_corners,
                                 feature_tracks, cameras, landmarks);

  return true;
}

// Localize a new camera in the map given a set of observed landmarks. We use
// pnp and ransac to localize the camera in the presence of outlier tracks.
// After finding an inlier set with pnp, we do non-linear refinement using all
// inliers and also update the set of inliers using the refined pose.
//
// shared_track_ids already contains those tracks which the new image shares
// with the landmarks (but some might be outliers).
//
// We return the refined pose and the set of track ids for all inliers.
//
// The inlier threshold is given in pixels. See also the opengv documentation on
// how to convert this to a ransac threshold:
// http://laurentkneip.github.io/opengv/page_how_to_use.html#sec_threshold
void localize_camera(
    const FrameCamId& fcid, const std::vector<TrackId>& shared_track_ids,
    const Calibration& calib_cam, const Corners& feature_corners,
    const FeatureTracks& feature_tracks, const Landmarks& landmarks,
    const double reprojection_error_pnp_inlier_threshold_pixel,
    Sophus::SE3d& T_w_c, std::vector<TrackId>& inlier_track_ids) {
  inlier_track_ids.clear();

  bearingVectors_t bvectors;
  points_t points;
  for (auto shared_tid : shared_track_ids) {
    FeatureTrack track = feature_tracks.at(shared_tid);
    auto f0 = track[fcid];

    auto c0 = feature_corners.at(fcid).corners[f0];

    bearingVector_t bvec(calib_cam.intrinsics[fcid.cam_id]->unproject(c0));
    point_t point = landmarks.at(shared_tid).p;
    bvectors.push_back(bvec);
    points.push_back(point);
  }

  absolute_pose::CentralAbsoluteAdapter adapter(bvectors, points);
  // TODO SHEET 4: Localize a new image in a given map
  sac::Ransac<sac_problems::absolute_pose::AbsolutePoseSacProblem> ransac;
  // create an AbsolutePoseSacProblem
  // (algorithm is selectable: KNEIP, GAO, or EPNP)
  std::shared_ptr<sac_problems::absolute_pose::AbsolutePoseSacProblem>
      absposeproblem_ptr(
          new sac_problems::absolute_pose::AbsolutePoseSacProblem(
              adapter,
              sac_problems::absolute_pose::AbsolutePoseSacProblem::KNEIP));
  // run ransac
  ransac.sac_model_ = absposeproblem_ptr;
  double f = 500;
  ransac.threshold_ =
      1.0 - cos(atan(reprojection_error_pnp_inlier_threshold_pixel / f));
  ransac.max_iterations_ = 5000;

  ransac.computeModel();
  if (ransac.inliers_.size() >= 3) {
    adapter.sett(ransac.model_coefficients_.block(0, 3, 3, 1));
    adapter.setR(ransac.model_coefficients_.block(0, 0, 3, 3));
    transformation_t refined_transformation =
        absolute_pose::optimize_nonlinear(adapter, ransac.inliers_);
    Eigen::Matrix4d pose = Eigen::Matrix4d::Identity();

    // refined_transformation.block(0,3,3,1) =  pose.block(0,3,3,1) ;
    ransac.sac_model_->selectWithinDistance(refined_transformation,
                                            ransac.threshold_, ransac.inliers_);

    pose.block(0, 0, 3, 3) = refined_transformation.block(0, 0, 3, 3);
    pose.block(0, 3, 3, 1) = refined_transformation.block(0, 3, 3, 1);
    T_w_c = Sophus::SE3d(pose);
    for (unsigned long i = 0; i < ransac.inliers_.size(); i++) {
      int ind = ransac.inliers_[i];
      inlier_track_ids.push_back(shared_track_ids[ind]);
    }
  }

  // get the result
}

struct BundleAdjustmentOptions {
  /// 0: silent, 1: ceres brief report (one line), 2: ceres full report
  int verbosity_level = 1;

  /// update intrinsics or keep fixed
  bool optimize_intrinsics = false;

  /// use huber robust norm or squared norm
  bool use_huber = true;

  /// parameter for huber loss (in pixel)
  double huber_parameter = 1.0;

  /// maximum number of solver iterations
  int max_num_iterations = 20;
};

// Run bundle adjustment to optimize cameras, points, and optionally intrinsics
void bundle_adjustment(const Corners& feature_corners,
                       const BundleAdjustmentOptions& options,
                       const std::set<FrameCamId>& fixed_cameras,
                       Calibration& calib_cam, Cameras& cameras,
                       Landmarks& landmarks) {
  ceres::Problem problem;
  ceres::HuberLoss* loss_function = NULL;
  // TODO SHEET 4: Setup optimization problem
  // ceres::HuberLoss *huber = new ceres::HuberLoss(options.huber_parameter);

  if (options.use_huber)
    loss_function = new ceres::HuberLoss(options.huber_parameter);

  for (auto& cam : cameras) {
    problem.AddParameterBlock(cam.second.T_w_c.data(),
                              Sophus::SE3d::num_parameters,
                              new Sophus::test::LocalParameterizationSE3);

    if (fixed_cameras.find(cam.first) != fixed_cameras.end()) {
      problem.SetParameterBlockConstant(cam.second.T_w_c.data());
    }
    double* intrParams =
        (double*)calib_cam.intrinsics[cam.first.cam_id]->getParam().data();
    problem.AddParameterBlock(intrParams, 8);
    if (options.optimize_intrinsics == false) {
      problem.SetParameterBlockConstant(intrParams);
    }
  }

  for (auto& lm : landmarks) {
    for (const auto& ob : lm.second.obs) {
      Eigen::Vector2d p_2d = feature_corners.at(ob.first).corners[ob.second];
      // std::cout << p_2d.transpose() << std::endl;
      ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<
          BundleAdjustmentReprojectionCostFunctor, 2,
          Sophus::SE3d::num_parameters, 3, 8>(
          new BundleAdjustmentReprojectionCostFunctor(
              p_2d, calib_cam.intrinsics[ob.first.cam_id]->name()));

      double* intrParams =
          (double*)calib_cam.intrinsics[ob.first.cam_id]->getParam().data();
      problem.AddResidualBlock(cost_function, loss_function,
                               cameras.at(ob.first).T_w_c.data(),
                               lm.second.p.data(), intrParams);
    }
  }

  // Solve
  // }

  ceres::Solver::Options ceres_options;
  ceres_options.max_num_iterations = options.max_num_iterations;
  ceres_options.linear_solver_type = ceres::SPARSE_SCHUR;
  ceres_options.num_threads = tbb::task_scheduler_init::default_num_threads();
  ceres::Solver::Summary summary;

  Solve(ceres_options, &problem, &summary);
  switch (options.verbosity_level) {
    // 0: silent
    case 1:
      std::cout << summary.BriefReport() << std::endl;
      break;
    case 2:
      std::cout << summary.FullReport() << std::endl;
      break;
  }
}

}  // namespace visnav
