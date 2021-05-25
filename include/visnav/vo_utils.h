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

#include <set>

#include <visnav/common_types.h>

#include <visnav/calibration.h>

#include <opengv/absolute_pose/CentralAbsoluteAdapter.hpp>
#include <opengv/absolute_pose/methods.hpp>
#include <opengv/relative_pose/CentralRelativeAdapter.hpp>
#include <opengv/sac/Ransac.hpp>
#include <opengv/sac_problems/absolute_pose/AbsolutePoseSacProblem.hpp>
#include <opengv/triangulation/methods.hpp>

namespace visnav {

void project_landmarks(
    const Sophus::SE3d& current_pose,
    const std::shared_ptr<AbstractCamera<double>>& cam,
    const Landmarks& landmarks, const double cam_z_threshold,
    std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>>&
        projected_points,
    std::vector<TrackId>& projected_track_ids) {
  projected_points.clear();
  projected_track_ids.clear();

  // TODO SHEET 5: project landmarks to the image plane using the current
  // locations of the cameras. Put 2d coordinates of the projected points into
  // projected_points and the corresponding id of the landmark into
  // projected_track_ids.
  for (auto& lm : landmarks) {
    auto point_in_cam = current_pose.inverse() * lm.second.p;
    if (point_in_cam.z() >= cam_z_threshold) {
      auto p_2d = cam->project(point_in_cam);
      if (p_2d.x() <= cam->width() && p_2d.y() <= cam->height() &&
          p_2d.x() >= 0 && p_2d.y() >= 0) {
        projected_points.push_back(p_2d);
        projected_track_ids.push_back(lm.first);
      }
    }
  }
}

int landmarkDistance(const Landmark& lm, const Corners& feature_corners,
                     const std::bitset<256>& descriptor) {
  int minDist = 256;

  for (const auto& ob : lm.obs) {
    std::bitset<256> desc =
        feature_corners.at(ob.first).corner_descriptors[ob.second];
    int dist = (desc ^ descriptor).count();
    if (dist < minDist) {
      minDist = dist;
    }
  }
  return minDist;
}

std::pair<int, int> arg_min2(int* arr, int size) {
  int idx =
      static_cast<int>(std::distance(arr, std::min_element(arr, arr + size)));
  int* minDist = std::min_element(arr, arr + size);
  std::pair<int, int> temp(idx, *minDist);
  return temp;
}
void find_matches_landmarks(
    const KeypointsData& kdl, const Landmarks& landmarks,
    const Corners& feature_corners,
    const std::vector<Eigen::Vector2d,
                      Eigen::aligned_allocator<Eigen::Vector2d>>&
        projected_points,
    const std::vector<TrackId>& projected_track_ids,
    const double match_max_dist_2d, const int feature_match_threshold,
    const double feature_match_dist_2_best, LandmarkMatchData& md) {
  md.matches.clear();
  // TODO SHEET 5: Find the matches between projected landmarks and detected
  // keypoints in the current frame. For every detected keypoint search for
  // matches inside a circle with radius match_max_dist_2d around the point
  // location. For every landmark the distance is the minimal distance between
  // the descriptor of the current point and descriptors of all observations of
  // the landmarks. The feature_match_threshold and feature_match_dist_2_best
  // should be used to filter outliers the same way as in exercise 3. You should
  // fill md.matches with <featureId,trackId> pairs for the successful matches
  // that pass all tests.
  int** distances = new int*[kdl.corners.size()];

  for (unsigned long int pt = 0; pt < kdl.corners.size(); pt++) {
    distances[pt] = new int[projected_points.size()];
    std::fill(distances[pt], distances[pt] + projected_points.size(), 0);
  }

  int** inv_distances = new int*[projected_points.size()];

  for (unsigned long int pt = 0; pt < kdl.corners.size(); pt++) {
    inv_distances[pt] = new int[kdl.corners.size()];
    std::fill(inv_distances[pt], inv_distances[pt] + kdl.corners.size(), 0);
  }

  for (unsigned long int i = 0; i < kdl.corners.size(); i++) {
    for (unsigned long int j = 0; j < projected_track_ids.size(); j++) {
      double dist = (projected_points[j] - kdl.corners[i]).norm();
      if (dist < match_max_dist_2d) {
        inv_distances[j][i] =
            landmarkDistance(landmarks.at(projected_track_ids[j]),
                             feature_corners, kdl.corner_descriptors[i]);
        distances[i][j] =
            landmarkDistance(landmarks.at(projected_track_ids[j]),
                             feature_corners, kdl.corner_descriptors[i]);

      } else {
        distances[i][j] = 256;
        inv_distances[j][i] = 256;
      }
    }
  }

  for (unsigned long int i = 0; i < kdl.corners.size(); i++) {
    std::pair<FeatureId, TrackId> m =
        arg_min2(distances[i], projected_track_ids.size());
    if (m.second >= feature_match_threshold) {
      continue;  // threshold check
    } else {
      int tempd = distances[i][m.first];
      distances[i][m.first] = 256;
      std::pair<FeatureId, TrackId> m3 =
          arg_min2(distances[i], projected_track_ids.size());
      distances[i][m.first] = tempd;
      // second best
      if (m3.second <= m.second * feature_match_dist_2_best) {
        continue;
      }  // else {
      md.matches.push_back(
          std::pair<FeatureId, TrackId>(i, projected_track_ids[m.first]));
      ;
    }
  }
}

// // for (long unsigned int i = 0; i < projected_track_ids.size(); i++) {
// //     delete distances[i];
// //   }

// //   for (long unsigned i = 0; i < kdl.corners.size(); i++) {
// //     delete inv_distances[i];
// //   }
// delete distances;
// delete inv_distances;

void localize_camera(const Sophus::SE3d& current_pose,
                     const std::shared_ptr<AbstractCamera<double>>& cam,
                     const KeypointsData& kdl, const Landmarks& landmarks,
                     const double reprojection_error_pnp_inlier_threshold_pixel,
                     LandmarkMatchData& md) {
  md.inliers.clear();

  // default to previous pose if not enough inliers
  md.T_w_c = current_pose;

  if (md.matches.size() < 4) {
    return;
  }

  // TODO SHEET 5: Find the pose (md.T_w_c) and the inliers (md.inliers) using
  // the landmark to keypoints matches and PnP. This should be similar to the
  // localize_camera in exercise 4 but in this exercise we don't explicitly have
  // tracks.

  bearingVectors_t bvectors;
  points_t points;
  for (auto match : md.matches) {
    auto c0 = kdl.corners[match.first];

    bearingVector_t bvec(cam->unproject(c0));
    point_t point = landmarks.at(match.second).p;
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
    md.T_w_c = Sophus::SE3d(pose);
    for (unsigned long i = 0; i < ransac.inliers_.size(); i++) {
      int ind = ransac.inliers_[i];
      md.inliers.push_back(md.matches[ind]);
    }
  }
}

void add_new_landmarks(const FrameCamId fcidl, const FrameCamId fcidr,
                       const KeypointsData& kdl, const KeypointsData& kdr,
                       const Calibration& calib_cam, const MatchData& md_stereo,
                       const LandmarkMatchData& md, Landmarks& landmarks,
                       TrackId& next_landmark_id) {
  // input should be stereo pair
  assert(fcidl.cam_id == 0);
  assert(fcidr.cam_id == 1);

  const Sophus::SE3d T_0_1 = calib_cam.T_i_c[0].inverse() * calib_cam.T_i_c[1];
  const Eigen::Vector3d t_0_1 = T_0_1.translation();
  const Eigen::Matrix3d R_0_1 = T_0_1.rotationMatrix();

  // TODO SHEET 5: Add new landmarks and observations. Here md_stereo contains
  // stereo matches for the current frame and md contains feature to landmark
  // matches for the left camera (camera 0). For all inlier feature to landmark
  // matches add the observations to the existing landmarks. If the left
  // camera's feature appears also in md_stereo.inliers, then add both
  // observations. For all inlier stereo observations that were not added to the
  // existing landmarks, triangulate and add new landmarks. Here
  // next_landmark_id is a running index of the landmarks, so after adding a new
  // landmark you should always increase next_landmark_id by 1.
  bearingVectors_t vec1;

  bearingVectors_t vec2;
  // run method 1
  std::vector<std::pair<int, int>> nonexistent_pairs;
  std::vector<TrackId> new_track_ids(md_stereo.inliers.size(), 1);
  for (auto& lm_match : md.inliers) {
    landmarks.at(lm_match.second)
        .obs.insert(std::make_pair(fcidl, lm_match.first));  // insert to left

    bool check = false;
    int idx = 0;

    for (auto& m : md_stereo.inliers) {
      if (kdl.corners[m.first] == kdl.corners[lm_match.first]) {
        check = true;
        new_track_ids[idx] = 0;
        break;
      }
      idx++;
    }
    if (check) {
      landmarks.at(lm_match.second)
          .obs.insert(std::make_pair(
              fcidr, md_stereo.inliers[idx].second));  // insert to right
    }
  }
  for (int i = 0; i < md_stereo.inliers.size(); i++) {
    if (new_track_ids[i] == 1)
      nonexistent_pairs.push_back(md_stereo.inliers[i]);
  }

  for (auto& st_match : nonexistent_pairs) {
    vec1.push_back(
        calib_cam.intrinsics[0]->unproject(kdl.corners[st_match.first]));
    vec2.push_back(
        calib_cam.intrinsics[1]->unproject(kdr.corners[st_match.second]));
  }

  relative_pose::CentralRelativeAdapter adapter(vec1, vec2, t_0_1, R_0_1);

  for (unsigned long j = 0; j < nonexistent_pairs.size(); j++) {
    point_t point = triangulation::triangulate(adapter, j);
    Landmark temp;
    temp.p = md.T_w_c * point;
    temp.obs.insert(std::make_pair(fcidl, nonexistent_pairs[j].first));
    temp.obs.insert(std::make_pair(fcidr, nonexistent_pairs[j].second));
    landmarks[next_landmark_id] = temp;
    next_landmark_id++;
  }
}

void remove_old_keyframes(const FrameCamId fcidl, const int max_num_kfs,
                          Cameras& cameras, Landmarks& landmarks,
                          Landmarks& old_landmarks,
                          std::set<FrameId>& kf_frames) {
  kf_frames.emplace(fcidl.frame_id);
  if (kf_frames.size() < max_num_kfs) {
    return;
  }
  // TODO SHEET 5: Remove old cameras and observations if the number of keyframe
  // pairs (left and right image is a pair) is larger than max_num_kfs. The ids
  // of all the keyframes that are currently in the optimization should be
  // stored in kf_frames. Removed keyframes should be removed from cameras and
  // landmarks with no left observations should be moved to old_landmarks.
  std::vector<int> kfs{kf_frames.begin(), kf_frames.end()};

  std::sort(kfs.begin(), kfs.end());
  int num_elems_to_delete = kfs.size() - max_num_kfs;
  std::vector<FrameCamId> deleted_cams;
  std::vector<TrackId> deleted_lms;
  for (int i = 0; i < num_elems_to_delete; i++) {
    auto elem = kf_frames.find(kfs[i]);
    if (elem != kf_frames.end()) kf_frames.erase(elem);
    for (auto& c : cameras) {
      if (c.first.frame_id == kfs[i]) {
        deleted_cams.push_back(c.first);
      }
    }
    for (auto& o : deleted_cams) {
      auto elem = cameras.find(o);
      if (elem != cameras.end()) cameras.erase(elem);
    }
    for (auto& lm : landmarks) {
      std::vector<FrameCamId> deleted_obs;

      for (auto& ob : lm.second.obs) {
        if (ob.first.frame_id == kfs[i]) {
          deleted_obs.push_back(ob.first);
        }
      }
      for (auto& o : deleted_obs) {
        auto elem = lm.second.obs.find(o);
        if (elem != lm.second.obs.end()) lm.second.obs.erase(elem);
      }
      if (lm.second.obs.size() - deleted_obs.size() == 0) {
        deleted_lms.push_back(lm.first);
        old_landmarks[lm.first] = lm.second;
      }
    }
    for (auto& o : deleted_lms) {
      auto elem = landmarks.find(o);
      if (elem != landmarks.end()) landmarks.erase(o);
    }
  }
}
}  // namespace visnav
