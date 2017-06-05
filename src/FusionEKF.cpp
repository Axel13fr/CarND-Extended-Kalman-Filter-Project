#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    H_laser_ << 1,0,0,0,
                0,1,0,0;
    Hj_ = MatrixXd(3, 4);
    MatrixXd Q = MatrixXd::Zero(4,4);

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
            0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
            0, 0.0009, 0,
            0, 0, 0.09;

    /**
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
    VectorXd x = VectorXd::Zero(4);
    MatrixXd mat4x4zero = MatrixXd::Zero(4,4);
    auto P = mat4x4zero;
    auto F = mat4x4zero;
    ekf_.Init(x,P,F,H_laser_,R_laser_,Q);

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::SetRadarMatrices()
{
    Hj_ = Tools::CalculateJacobian(ekf_.x_);
    ekf_.setH(Hj_);
    ekf_.setR(R_radar_);
}

void FusionEKF::SetLaserMatrices()
{
    ekf_.setH(H_laser_);
    ekf_.setR(R_laser_);
}

void FusionEKF::UpdateFandQ(const MeasurementPackage &measurement_pack)
{
    constexpr double MICROSECS_TO_SECS = 1/1000000.0;
    auto delta_t = (measurement_pack.timestamp_us_ - previous_timestamp_)*MICROSECS_TO_SECS;
    // F update
    ekf_.F_ << 1,0,delta_t,0,
            0,1,0,delta_t,
            0,0,1,0,
            0,0,0,1;

    // Q update
    auto delta4 = std::pow(delta_t,4) / 4.0;
    auto delta3 = std::pow(delta_t,3) / 2.0;
    auto delta2 = delta_t*delta_t;
    constexpr auto NOISE_AX = 9.0;
    constexpr auto NOISE_AY = 9.0;
    ekf_.Q_ <<  delta4*NOISE_AX ,      0   ,    delta3*NOISE_AX,      0,
                 0      ,       delta4*NOISE_AY,    0     ,     delta3*NOISE_AY,
                delta3*NOISE_AX,       0 ,     delta2*NOISE_AX,    0,
                  0,           delta3*NOISE_AY,      0 ,        delta2*NOISE_AY;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


    /*****************************************************************************
   *  Initialization
   ****************************************************************************/

    if (!is_initialized_) {
        /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
        // first measurement
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;
        ekf_.P_ = MatrixXd(4, 4);
        // Initial cov set to 2m for position and 1m/s for speed
        cout << "Init P cov matrix:" << endl;
        ekf_.P_ << 3,0,0,0,
                0,3,0,0,
                0,0,3,0,
                0,0,0,3;
        cout << ekf_.P_ << endl;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
            ekf_.x_ = Tools::CalculatePosFromRadar(measurement_pack.raw_measurements_);
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
      Initialize state.
      */
            ekf_.x_ <<  measurement_pack.raw_measurements_(0),
                    measurement_pack.raw_measurements_(1),
                    0.1,0.1;
        }

        previous_timestamp_ = measurement_pack.timestamp_us_;
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
   *  Prediction
   ****************************************************************************/

    /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
    UpdateFandQ(measurement_pack);
    ekf_.Predict();

    /*****************************************************************************
   *  Update
   ****************************************************************************/

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        SetRadarMatrices();
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
        previous_timestamp_ = measurement_pack.timestamp_us_;

    } else {
        // Laser updates
        SetLaserMatrices();
        ekf_.Update(measurement_pack.raw_measurements_);
        previous_timestamp_ = measurement_pack.timestamp_us_;
    }

    // print the output
    cout << "Fused pos x_ = " << ekf_.x_(0) << " y_ = " << ekf_.x_(1)  << " Vx_ = " << ekf_.x_(2) << " Vy_ = " << ekf_.x_(3) << endl;
    cout << "sig_x = " << ekf_.P_(0,0) << " sig_y = " << ekf_.P_(1,1) <<  endl;
}
