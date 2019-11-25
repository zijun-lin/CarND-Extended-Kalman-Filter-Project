#include "kalman_filter.h"
#include <iostream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
    VectorXd y = z - H_ * x_;
    MatrixXd S = H_ * P_ * H_.transpose() + R_;
    MatrixXd K = P_ * H_.transpose() * S.inverse();
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);

    x_ = x_ + K * y;
    P_ = (I - K *H_)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
    // state parameters
    double px = x_(0);
    double py = x_(1);
    double vx = x_(2);
    double vy = x_(3);
    // pre-compute a set of terms to avoid repeated calculation
    float rho = sqrtf(px*px + py*py);
    float phi = atan2(py, px);
    float rho_dot;
    if(fabs(rho) < 0.000001){
        rho_dot = 0.0;
    } else{
        rho_dot =(px * vx + py * vy)/rho;
    }
    VectorXd hx(3);
    hx << rho, phi, rho_dot;


    VectorXd y = z - hx;
    float angle = y(1);
    const float PI = 3.14159265;
    if (angle > PI || angle < -PI) {
        float new_angle = std::fmod(angle, 2.0 * PI);
        angle = new_angle < 0 ? new_angle + 2.0 * PI : new_angle;
        y(1) = angle;
        std::cout << " Normalizing Angles " << std::endl;
    }

    MatrixXd S = H_ * P_ * H_.transpose() + R_;
    MatrixXd K = P_ * H_.transpose() * S.inverse();
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);

    x_ = x_ + K * y;
    P_ = (I - K * H_) * P_;
}

// behavioral cloning