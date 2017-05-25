#include <iostream>
#include <assert.h>
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

bool Tools::CalculateJacobian(const VectorXd& x_state,MatrixXd &Hj) {
  /**
    * Calculate a Jacobian here.
  */
  assert(Hj.rows()==3 && Hj.cols()==4);
  assert(x_state.size() ==4);
	// initialize the Jacobian Matrix

	//recover state parameters
	double px = x_state(0);
	double py = x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	double p2= px * px + py * py; //square of p
	double p = sqrt(p2);           //length of vector or distance
	double p3 = p * p2;           //third power of p

	//check division by zero
	if (p == 0) {
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return false;
	}

	//compute the Jacobian matrix
	Hj << px / p, 		                  py / p, 		                   0, 		 0,
			  -py / p2,                     px / p2,                       0,      0,
			  py* (vx * py - vy * px) / p3, px * (px * vy - py * vx) / p3, px / p, py/ p;

	return true;
}

double Tools::NormalizeMinMax(double x, double min, double max){

  assert(min <max);
  double width = max - min;
  double offset = x - min;
  return offset - std::floor(offset/width)*width + min;
}
