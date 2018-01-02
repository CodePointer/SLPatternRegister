//
// Created by pointer on 17-12-15.
//

#ifndef PATTERNREFINEMENT_INTENSITY_FUNCTOR_H
#define PATTERNREFINEMENT_INTENSITY_FUNCTOR_H

#include <ceres/ceres.h>

class IntensityFunctor {
public:
  IntensityFunctor(Eigen::Matrix<double, 3, Eigen::Dynamic> point_set,
                   double img_k);

  template <class T>
  bool operator()(const T* const C_p, const T* const Sigma,
                  T* residuals) const {
    long point_num = point_set_.cols();
    T img_k_head = T(0);
    for(int i = 0; i < point_num; i++) {
      int color_num = int(this->point_set_(0, i));
      double d_x = this->point_set_(1, i);
      double d_y = this->point_set_(2, i);
      T exp_term = ceres::exp(-T(d_x*d_x)/(T(2)*Sigma[0]*Sigma[0])-T(d_y*d_y)/(T(2)*Sigma[0]*Sigma[0]));
      img_k_head += T(1) / (T(2*M_PI) * Sigma[0] * Sigma[0]) * exp_term * C_p[color_num];
    }
    residuals[0] = img_k_ - img_k_head;
    return true;
  }


  // [[color_1; d_x_1; d_y_1], [color_2; d_x_1; d_y_1], ...]
  Eigen::Matrix<double, 3, Eigen::Dynamic> point_set_;
  double img_k_;
};


#endif //PATTERNREFINEMENT_INTENSITY_FUNCTOR_H
