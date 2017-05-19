#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE here.
  */

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if(estimations.size() == 0 || estimations.size() != ground_truth.size()){
	  cout << "Invalid estimation or ground_truth data" << endl;
	  return rmse;
  }

  unsigned int size = estimations.size();

  for(unsigned int i=0; i<size; i++){
	  VectorXd residual = estimations[i] - ground_truth[i];
	  residual = residual.array() * residual.array();
	  rmse += residual;
  }
  //calculate the mean
  rmse = rmse/size;

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian here.
  */

	MatrixXd Hj(3,4);
	// initialize the Jacobian Matrix
	Hj << 0,0,0,0,
			  0,0,0,0,
				0,0,0,0;
	//recover state parameters
	double px = x_state(0);
	double py = x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	double p2= px * px + py * py; //square of p
	double p = sqrt(p);           //length of vector or distance
	double p3 = p * p2;           //third power of p

	//check division by zero
	if (p2 == 0) {
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj;
	}

	//compute the Jacobian matrix
	Hj << px / p, 		                  py / p, 		                   0, 		 0,
			  -py / p2,                     px / p2,                       0,      0,
			  py* (vx * py - vy * px) / p3, px * (px * vy - py * vx) / p3, px / p, py/ p;

	return Hj;
}
