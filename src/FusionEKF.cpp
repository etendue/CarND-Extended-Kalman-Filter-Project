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
  Hj_ = MatrixXd(3, 4);

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
  MatrixXd P = MatrixXd(4,4);
  MatrixXd F = MatrixXd(4,4);
  MatrixXd Q = MatrixXd(4,4);

	P << 1, 0, 0, 0,
			 0, 1, 0, 0,
			 0, 0, 1000, 0,
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
	ekf_.Init(ekf_.x_,P,F,ekf_.H_,ekf_.R_,Q);

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
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
			ekf_.H_ = Hj_;
			ekf_.R_ = R_radar_;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
    	if (measurement_pack.raw_measurements_.size()!=4){
    		cout << "Invalid measurement data !!!" << endl;
    		return ;
    	}
    	//initialize x
    	x << measurement_pack.raw_measurements_[0],
    			 measurement_pack.raw_measurements_[1],
					 measurement_pack.raw_measurements_[2],
					 measurement_pack.raw_measurements_[3];

    	//assign the measurement matrix and measurement noise
    	ekf_.x_ = x;
    	ekf_.H_ = H_laser_;
    	ekf_.R_ = R_laser_;

    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  double dt = (measurement_pack.timestamp_ - previous_timestamp_)/1.0e6;
  double dt2 = dt*dt; // delta t square
  double dt3 = dt*dt2;// delta t power 3
  double dt4 = dt2 * dt2; // delta t power 4
  //set current timestamp as last timestamp for next cycle
  previous_timestamp_ = measurement_pack.timestamp_;
  //copy the last value
  MatrixXd F = ekf_.F_;
  MatrixXd Q = ekf_.Q_;
  // update prediction F matrix with dt value
  F(0,2) = F(1,3) = dt;
  // update process noise matrix with dt value
  Q <<  dt4 * 9/4, 0, dt3 * 9/2, 0,
        0, dt4 * 9/4, 0, dt3 * 9/2,
				dt3 * 9/2, 0, dt2 * 9, 0,
        0, dt3 * 9/2, 0, dt2 * 9;
  // assign the updated F and Q to kalman filter
  ekf_.F_ = F;
  ekf_.Q_ = Q;
  ekf_.Predict();

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

		VectorXd z = VectorXd(3);
		z << measurement_pack.raw_measurements_[0],
				 measurement_pack.raw_measurements_[1],
				 measurement_pack.raw_measurements_[2];

		// get Jacob Matrix
		ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
		ekf_.UpdateEKF(z);
  } else {
    // Laser updates
		if (measurement_pack.raw_measurements_.size() != 4) {
			cout << "Invalid measurement data !!!" << endl;
			return;
		}
		//initialize x
		VectorXd z = VectorXd(4);
		z << measurement_pack.raw_measurements_[0],
				 measurement_pack.raw_measurements_[1],
				 measurement_pack.raw_measurements_[2],
				 measurement_pack.raw_measurements_[3];

		//assign the measurement matrix and measurement noise
		ekf_.Update(z);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
