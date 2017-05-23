#include "FusionEKF.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

static double noise_ax = 9.0;
static double noise_ay = 9.0;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;
  Qv_ = VectorXd(2);
  Qv_<< noise_ax,noise_ay;

  // initializing matrices
  MatrixXd R_laser_ = MatrixXd(2, 2);
  MatrixXd R_radar_ = MatrixXd(3, 3);
  MatrixXd H_laser_ = MatrixXd(2, 4);
  MatrixXd Hj_ = MatrixXd(3, 4);


  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   TODO:
   * Finish initializing the FusionEKF.
   * Set the process and measurement noises
   */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  MatrixXd P = MatrixXd(4, 4);
  MatrixXd F = MatrixXd(4, 4);
  MatrixXd Q = MatrixXd(4, 4);

  P << 1, 0, 0, 0,
       0, 1, 0, 0,
       0, 0, 1000,0,
       0, 0, 0, 1000;

  F << 1, 0, 1, 0,
       0, 1, 0, 1,
       0, 0, 1, 0,
       0, 0, 0, 1;

  Q << 0, 0, 0, 0,
       0, 0, 0, 0,
       0, 0, 0, 0,
       0, 0, 0, 0;



  // only initialize the P, and F, others are not known now
  ekf_.Init(ekf_.x_, P, F, H_laser_,Hj_,R_laser_,R_radar_, Q);

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {
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
    VectorXd x = VectorXd(4);
    x << 0, 0, 0, 0;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
       Convert radar from polar to cartesian coordinates and initialize state.
       */
      if (measurement_pack.raw_measurements_.size() != 3) {
        cout << "Invalid measurement data !!!" << endl;
        return;
      }

      double range = measurement_pack.raw_measurements_[0];
      double bearing = measurement_pack.raw_measurements_[1];
      double range_rate = measurement_pack.raw_measurements_[2];
      x[0] = range * cos(bearing);
      x[1] = range * sin(bearing);
      x[2] = range_rate * cos(bearing);
      x[3] = range_rate * sin(bearing);

      ekf_.x_ = x;
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
       Initialize state.
       */
      if (measurement_pack.raw_measurements_.size() != 2) {
        cout << "Invalid measurement data !!!" << endl;
        return;
      }
      //initialize x
      x << measurement_pack.raw_measurements_[0], measurement_pack
          .raw_measurements_[1], 0, 0;
      //assign the measurement matrix
      ekf_.x_ = x;
    }

    //set current timestamp as last timestamp for next cycle
    previous_timestamp_ = measurement_pack.timestamp_;

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
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1.0e6;
  //set current timestamp as last timestamp for next cycle
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.Predict(dt,Qv_);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
   * Use the sensor type to perform the update step.
   * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    if (measurement_pack.raw_measurements_.size() != 3) {
      cout << "Invalid measurement data !!!" << endl;
      return;
    }

    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
    if (measurement_pack.raw_measurements_.size() != 2) {
      cout << "Invalid measurement data !!!" << endl;
      return;
    }
    //assign the measurement matrix and measurement noise
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
