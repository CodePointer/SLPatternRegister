#include <iostream>
#include <opencv2/opencv.hpp>
#include <Eigen/Eigen>
#include <ceres/ceres.h>
#include <vector>
#include "intensity_functor.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ImgMatrix;

const int kHeight = 1024;
const int kWidth = 1280;
const int kProHeight = 800;
const int kProWidth = 1280;

cv::Mat LoadTxtFileToMat(std::string file_name, int kHeight, int kWidth) {
  cv::Mat result(kHeight, kWidth, CV_64FC1);
  std::fstream file(file_name, std::ios::in);
  for (int h = 0; h < kHeight; h++) {
    for (int w = 0; w < kWidth; w++) {
      double tmp;
      file >> tmp;
      result.at<double>(h, w) = tmp;
    }
  }
  file.close();
  return result;
}

int main() {
  // Read part
  cv::Mat img = cv::imread("img.png", cv::IMREAD_GRAYSCALE);
  cv::Mat pattern = cv::imread("pattern.png", cv::IMREAD_GRAYSCALE);
  cv::Mat x_pro_mat = LoadTxtFileToMat("xpro_mat0.txt", kHeight, kWidth);
  cv::Mat y_pro_mat = LoadTxtFileToMat("ypro_mat0.txt", kHeight, kWidth);
  cv::Mat mask_mat = LoadTxtFileToMat("mask.txt", kHeight, kWidth);

  // ceres initialization
  double color_set[2] = {60, 628};
  double sigma = 1.0;
  ceres::Problem * problem = new ceres::Problem;
  ceres::Solver::Options options;
  options.gradient_tolerance = 1e-8;
  options.function_tolerance = 1e-8;
  options.linear_solver_type = ceres::DENSE_QR;
  options.max_num_iterations = 70;
  options.minimizer_progress_to_stdout = true;
  ceres::Solver::Summary summary;

  // Every pixel
  int kSearchRad = 10;
  int block_num = 0;
  for (int h = 0; h < kHeight; h++) {
    for (int w = 0; w < kWidth; w++) {
      double mask_valid = mask_mat.at<double>(h, w);
      if (mask_valid < 0.5) {
        continue;
      }
      double x_pro = x_pro_mat.at<double>(h, w) - 5;
      double y_pro = y_pro_mat.at<double>(h, w);
      if ((x_pro < 0) || (y_pro < 0)) {
        continue;
      }
      std::vector<cv::Point2d> vec_h_w;
      for (int d_h = -kSearchRad; d_h <= kSearchRad; d_h++) {
        for (int d_w = -kSearchRad; d_w <= kSearchRad; d_w++) {
          int h_n = int(y_pro) + d_h;
          int w_n = int(x_pro) + d_w;
          if ((h_n >= 0) && (h_n < kProHeight)
              && (w_n >= 0) && (w_n < kProWidth)) {
            if ((h_n % 2 == 0) && (w_n % 2 == 0)) {
              vec_h_w.emplace_back(cv::Point2d(w_n + 0.5, h_n + 0.5));
            }
          }
        }
      }
      Eigen::Matrix<double, 3, Eigen::Dynamic> point_set(3, vec_h_w.size());
      for (int i = 0; i < vec_h_w.size(); i++) {
        double x_grid = vec_h_w[i].x;
        double y_grid = vec_h_w[i].y;
        uchar tmp_color = pattern.at<uchar>(int(y_grid-0.5), int(x_grid-0.5));
        if (tmp_color == 255) {
          point_set(0, i) = 1;
        } else {
          point_set(0, i) = 0;
        }
        point_set(1, i) = std::abs(x_grid - x_pro);
        point_set(2, i) = std::abs(y_grid - y_pro);
      }
      double img_k = double(img.at<uchar>(h, w));

      ceres::CostFunction * cost_fun =
          new ceres::AutoDiffCostFunction<IntensityFunctor, 1, 2, 1>(
              new IntensityFunctor(point_set, img_k));
      problem->AddResidualBlock(cost_fun, nullptr, color_set, &sigma);
      block_num++;
    }
  }
  std::cout << block_num << std::endl;

  // Solve
  ceres::Solve(options, problem, &summary);
  std::cout << summary.FullReport() << std::endl;

  std::cout << color_set[0] << std::endl;
  std::cout << color_set[1] << std::endl;
  std::cout << sigma << std::endl;

  // Generated image
  cv::Mat img_obs(kHeight, kWidth, CV_8UC1);
  img_obs.setTo(0);
  for (int h = 0; h < kHeight; h++) {
    for (int w = 0; w < kWidth; w++) {
      double mask_valid = mask_mat.at<double>(h, w);
      if (mask_valid < 0.5) {
        continue;
      }
      double x_pro = x_pro_mat.at<double>(h, w) - 5;
      double y_pro = y_pro_mat.at<double>(h, w);
      if ((x_pro < 0) || (y_pro < 0)) {
        continue;
      }
      double intensity = 0.0;
      for (int d_h = -kSearchRad; d_h <= kSearchRad; d_h++) {
        for (int d_w = -kSearchRad; d_w <= kSearchRad; d_w++) {
          int h_n = int(y_pro) + d_h;
          int w_n = int(x_pro) + d_w;
          if ((h_n >= 0) && (h_n < kProHeight)
              && (w_n >= 0) && (w_n < kProWidth)) {
            if ((h_n % 2 == 0) && (w_n % 2 == 0)) {
              double color;
              if (pattern.at<uchar>(h_n, w_n) == 255) {
                color = color_set[1];
              } else {
                color = color_set[0];
              }
              double dx = std::abs(w_n + 0.5 - x_pro);
              double dy = std::abs(h_n + 0.5 - y_pro);
              double exp_term = ceres::exp(-(dx*dx)/(2*sigma*sigma)-(dy*dy)/(2*sigma*sigma));
              intensity += color / (2*M_PI*sigma*sigma) * exp_term;
            }
          }
        }
      }
      img_obs.at<uchar>(h, w) = uchar(intensity);
    }
  }
  cv::imwrite("img_obs.png", img_obs);

  return 0;
}