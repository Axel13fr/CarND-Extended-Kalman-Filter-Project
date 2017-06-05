#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

void KalmanFilter::setH(const Eigen::MatrixXd &H)
{
    H_ = H;
}

void KalmanFilter::setR(const Eigen::MatrixXd &R)
{
    R_ = R;
}

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in; // Transforms from state space to measurement space (depends on meas. type !)
    R_ = R_in; // measurment cov matrix (depends on meas. type !)
    Q_ = Q_in;
}

void KalmanFilter::Predict() {
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

    /*
     * KF Measurement update step : laser case
     * H and R must be set appropriately before the call !
     */
    VectorXd y = z - H_ * x_;
    MatrixXd S = H_ * P_ * H_.transpose() + R_;
    MatrixXd K =  P_ * H_.transpose() * S.inverse();

    //new state
    x_ = x_ + (K * y);
    auto I = MatrixXd::Identity(4, 4);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    /*
     * KF Measurement update step : radar case
     * H and R must be set appropriately before the call !
     */
    VectorXd y = z - H_ * x_;
    MatrixXd S = H_ * P_ * H_.transpose() + R_;
    auto K =  P_ * H_.transpose() * S.inverse();

    //new state
    x_ = x_ + (K * y);
    auto I = MatrixXd::Identity(3, 3);
    P_ = (I - K * H_) * P_;
}
