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

#include <sophus/se3.hpp>

#include <visnav/common_types.h>

namespace visnav {
// Taylor approximations from wolphram alpha
// sin(x)/x
template <class T>
T approximateSinTheta(T x) {
  T x2 = x * x;
  T x4 = x2 * x2;
  return 1.0f - (x2) / 6 + x4 / 120 + x4 * x2;
}

template <class T>
T approximateCosTheta(T x) {
  T x2 = x * x;
  T x4 = x2 * x2;
  return 0.5f - (x2) / 24 + (x2 * x2) / 720 + x4 * x2;
}

template <class T>
T approximateCosThetaSq(T x) {
  T x2 = x * x;
  T x4 = x2 * x2;

  return 0.5f - x2 / 24 + x4 / 720 + x4 * x2;
}
template <class T>
T approximateSinThetaCube(T x) {
  T x2 = x * x;
  T x4 = x2 * x2;

  return 1.0f / 6 - x2 / 120 + x4 / 5040 + x4 * x2;
}

template <class T>
T approximateSinCosThetaSq(T x) {
  T x2 = x * x;
  T x4 = x2 * x2;

  return 1.0f / 12 + x2 / 720 + x4 / 30240 + x4 * x2;
}
// Implement exp for SO(3)
template <class T>
Eigen::Matrix<T, 3, 3> user_implemented_expmap(
    const Eigen::Matrix<T, 3, 1>& xi) {
  // TODO SHEET 1: implement
  Eigen::Matrix<T, 3, 3> w_hat;
  w_hat << 0, -xi(2), xi(1), xi(2), 0, -xi(0), -xi(1), xi(0), 0;

  T theta = xi.norm();
  Eigen::Matrix<T, 3, 3> res;
  if (abs(theta) < 1e-6) {
    // Eigen::Matrix<T, 3, 3>::Identity()
    res = Eigen::Matrix<T, 3, 3>::Identity() +
          w_hat * (approximateSinTheta(theta)) +
          (approximateCosTheta(theta)) * w_hat * w_hat;
  } else {
    res = Eigen::Matrix<T, 3, 3>::Identity() + w_hat * sin(theta) / theta +
          ((1 - cos(theta)) / (theta * theta)) * w_hat * w_hat;
  }
  return res;
}

// Implement log for SO(3)
template <class T>
Eigen::Matrix<T, 3, 1> user_implemented_logmap(
    const Eigen::Matrix<T, 3, 3>& mat) {
  // TODO SHEET 1: implement
  T theta = acos((mat.trace() - 1) / 2);

  Eigen::Matrix<T, 3, 1> w;
  w << mat(2, 1) - mat(1, 2), mat(0, 2) - mat(2, 0), mat(1, 0) - mat(0, 1);

  if (abs(theta) < 1e-6) {
    w = w * approximateSinTheta(theta) * 1 /
        2.0f;  // Eigen::Matrix<T, 3, 1>::Zero();
  } else {
    w = (theta / (2 * sin(theta))) * w;
  }
  return w;
}

// Implement exp for SE(3)
template <class T>
Eigen::Matrix<T, 4, 4> user_implemented_expmap(
    const Eigen::Matrix<T, 6, 1>& xi) {
  // TODO SHEET 1: implement

  Eigen::Matrix<T, 3, 1> w;
  w << xi(3), xi(4), xi(5);

  Eigen::Matrix<T, 3, 1> v;
  v << xi(0), xi(1), xi(2);
  T theta = w.norm();
  Eigen::Matrix<T, 3, 3> w_hat;
  w_hat << 0, -w(2), w(1), w(2), 0, -w(0), -w(1), w(0), 0;

  Eigen::Matrix<T, 3, 3> rot = user_implemented_expmap(w);
  Eigen::Matrix<T, 3, 3> J;
  if (abs(theta) < 1e-6) {
    J = Eigen::Matrix<T, 3, 3>::Identity() +
        (approximateCosThetaSq(theta)) * w_hat +
        (approximateSinThetaCube(theta)) * w_hat * w_hat;
  } else {
    J = Eigen::Matrix<T, 3, 3>::Identity() +
        ((1 - cos(theta)) / (theta * theta)) * w_hat +
        ((theta - sin(theta)) / (theta * theta * theta)) * w_hat * w_hat;
  }
  Eigen::Matrix<T, 4, 4> res = Eigen::Matrix<T, 4, 4>::Identity();
  res.block(0, 0, 3, 3) = rot;
  res.block(0, 3, 3, 1) = J * v;

  return res;
}

// Implement log for SE(3)
template <class T>
Eigen::Matrix<T, 6, 1> user_implemented_logmap(
    const Eigen::Matrix<T, 4, 4>& mat) {
  // TODO SHEET 1: implement
  Eigen::Matrix<T, 3, 3> rot = mat.block(0, 0, 3, 3);
  Eigen::Matrix<T, 3, 1> t = mat.block(0, 3, 3, 1);
  Eigen::Matrix<T, 3, 1> w = user_implemented_logmap(rot);
  T theta = w.norm();
  Eigen::Matrix<T, 3, 3> w_hat;
  w_hat << 0, -w(2), w(1), w(2), 0, -w(0), -w(1), w(0), 0;
  Eigen::Matrix<T, 3, 3> J_inv;
  if (abs(theta) < 1e-6) {
    J_inv = Eigen::Matrix<T, 3, 3>::Identity() - (1.0f / 2.0f) * w_hat +
            (approximateSinCosThetaSq(theta)) * w_hat * w_hat;

  } else {
    J_inv = Eigen::Matrix<T, 3, 3>::Identity() - (1.0f / 2.0f) * w_hat +
            (1.0f / (theta * theta) -
             ((1.0f + cos(theta)) / (2.0f * theta * sin(theta)))) *
                w_hat * w_hat;
  }
  Eigen::Matrix<T, 3, 1> v = J_inv * t;
  // J_inv = Eigen::Matrix<T,3,3>::Identity();
  Eigen::Matrix<T, 6, 1> res = Eigen::Matrix<T, 6, 1>::Zero();
  res.block(0, 0, 3, 1) = v;
  res.block(3, 0, 3, 1) = w;
  return res;
}

}  // namespace visnav
