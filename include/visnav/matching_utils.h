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

#include <bitset>
#include <set>

#include <Eigen/Dense>
#include <sophus/se3.hpp>

#include <opengv/relative_pose/CentralRelativeAdapter.hpp>
#include <opengv/relative_pose/methods.hpp>
#include <opengv/sac/Ransac.hpp>
#include <opengv/sac_problems/relative_pose/CentralRelativePoseSacProblem.hpp>

#include <visnav/camera_models.h>
#include <visnav/common_types.h>
using namespace opengv;
namespace visnav {

void computeEssential(const Sophus::SE3d& T_0_1, Eigen::Matrix3d& E) {
  const Eigen::Vector3d t_0_1 = T_0_1.translation();
  const Eigen::Matrix3d R_0_1 = T_0_1.rotationMatrix();
  Eigen::Vector3d t = t_0_1.normalized();
  // TODO SHEET 3: compute essential matrix
  Eigen::Matrix3d t_hat = Sophus::SO3d::hat(t);
  E = t_hat * R_0_1;
  // std::cout << E << std::endl;
}

void findInliersEssential(const KeypointsData& kd1, const KeypointsData& kd2,
                          const std::shared_ptr<AbstractCamera<double>>& cam1,
                          const std::shared_ptr<AbstractCamera<double>>& cam2,
                          Eigen::Matrix3d& E, double epipolar_error_threshold,
                          MatchData& md) {
  md.inliers.clear();
  for (size_t j = 0; j < md.matches.size(); j++) {
    const Eigen::Vector2d p0_2d = kd1.corners[md.matches[j].first];
    const Eigen::Vector2d p1_2d = kd2.corners[md.matches[j].second];

    Eigen::Vector3d x1_homo = cam1->unproject(p0_2d);
    Eigen::Vector3d x2_homo = cam2->unproject(p1_2d);

    double epc = abs(x1_homo.transpose() * E * x2_homo);

    if (epc <= epipolar_error_threshold) {
      md.inliers.push_back(md.matches[j]);
    }
  }
}

void findInliersRansac(const KeypointsData& kd1, const KeypointsData& kd2,
                       const std::shared_ptr<AbstractCamera<double>>& cam1,
                       const std::shared_ptr<AbstractCamera<double>>& cam2,
                       const double ransac_thresh, const int ransac_min_inliers,
                       MatchData& md) {
  md.inliers.clear();
  md.T_i_j = Sophus::SE3d();
  bearingVectors_t vec1;
  bearingVectors_t vec2;
  for (size_t j = 0; j < md.matches.size(); j++) {
    const Eigen::Vector2d p0_2d = kd1.corners[md.matches[j].first];
    const Eigen::Vector2d p1_2d = kd2.corners[md.matches[j].second];

    Eigen::Vector3d x1_homo = cam1->unproject(p0_2d);
    Eigen::Vector3d x2_homo = cam2->unproject(p1_2d);
    vec1.push_back(x1_homo);
    vec2.push_back(x2_homo);
  }

  relative_pose::CentralRelativeAdapter adapter(vec1, vec2);
  // create a RANSAC object
  sac::Ransac<sac_problems::relative_pose::CentralRelativePoseSacProblem>
      ransac;
  // create a CentralRelativePoseSacProblem
  // (set algorithm to STEWENIUS, NISTER, SEVENPT, or EIGHTPT)
  std::shared_ptr<sac_problems::relative_pose::CentralRelativePoseSacProblem>
      relposeproblem_ptr(
          new sac_problems::relative_pose::CentralRelativePoseSacProblem(
              adapter, sac_problems::relative_pose::
                           CentralRelativePoseSacProblem::NISTER));
  // run ransac
  ransac.sac_model_ = relposeproblem_ptr;
  ransac.threshold_ = ransac_thresh;
  ransac.max_iterations_ = 5000;
  ransac.computeModel();
  // get the result
  if (ransac.inliers_.size() >= ransac_min_inliers) {
    transformation_t refined_transformation =
        relative_pose::optimize_nonlinear(adapter, ransac.inliers_);
    Eigen::Matrix4d pose = Eigen::Matrix4d::Identity();

    // refined_transformation.block(0,3,3,1) =  pose.block(0,3,3,1) ;
    ransac.sac_model_->selectWithinDistance(refined_transformation,
                                            ransac_thresh, ransac.inliers_);
    pose.block(0, 0, 3, 3) = refined_transformation.block(0, 0, 3, 3);
    pose.block(0, 3, 3, 1) =
        refined_transformation.block(0, 3, 3, 1).normalized();
    md.T_i_j = Sophus::SE3d(pose);
    for (int i = 0; i < ransac.inliers_.size(); i++) {
      int ind = ransac.inliers_[i];
      md.inliers.push_back(md.matches[ind]);
    }
  }
}

}  // namespace visnav
