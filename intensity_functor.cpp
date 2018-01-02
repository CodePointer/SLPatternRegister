//
// Created by pointer on 17-12-15.
//

#include "intensity_functor.h"

#include <utility>

IntensityFunctor::IntensityFunctor(
    Eigen::Matrix<double, 3, Eigen::Dynamic> point_set, double img_k) {
  this->point_set_ = std::move(point_set);
  this->img_k_ = img_k;
}
