#include "kalman_filter.h"
#include <iostream>
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

constexpr double pi() { return acos(-1.0);}

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
      MatrixXd &H_in,MatrixXd &Hj_in,MatrixXd &Rl_in,MatrixXd &Rr_in, MatrixXd &Q_in)
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  Hj_ = Hj_in;
  Rl_ = Rl_in;
  Rr_ = Rr_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict(const double dt, const VectorXd &noise) {
  /**
    * predict the state
  */
  // update prediction F matrix with dt value
  F_(0, 2) = F_(1, 3) = dt;

  double dt2 = dt * dt;  // delta t square
  double dt3 = dt * dt2;  // delta t power 3
  double dt4 = dt2 * dt2;  // delta t power 4

  if(noise.size()!=2){
    cout << "Invalid Process noise vector!!!" <<endl;
    return;
  }
  double noise_ax = noise[0];
  double noise_ay = noise[1];

    // update process noise matrix with dt value
  Q_ << dt4 * noise_ax / 4, 0, dt3 * noise_ax / 2, 0,
        0, dt4 * noise_ay/ 4, 0, dt3 * noise_ay / 2,
        dt3 * noise_ax / 2, 0, dt2 * noise_ax,
        0, 0, dt3 * noise_ay / 2, 0, dt2 * noise_ay;
	//x' = F x + u
	x_ = F_*x_;
	//P' = F P Ft + Q
	P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */
	VectorXd y = z - H_ * x_;
	MatrixXd S = H_ * P_ * H_.transpose() + Rl_;
	MatrixXd K =  P_ * H_.transpose() * S.inverse();

	//new estimate
	x_ = x_ + (K * y);
	unsigned int x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */
  double px = x_[0];
  double py = x_[1];
  double vx = x_[2];
  double vy = x_[3];

  // calculate h(x')
  double hx_range =  sqrt(px * px +py*py);
  if(hx_range == 0)
  {
    cout << "predict range is zero" << endl;
    return;
  }
  double hx_bearing = atan2(py,px);
  double hx_range_rate = (px*vx + py*vy)/hx_range;

  // calcuate y = z - h(x')
  VectorXd hx = VectorXd(3);
  hx << hx_range,hx_bearing,hx_range_rate;

  double mea_bearing = tools.NormalizeMinMax(z[1],-pi(),pi());
  cout <<" pre:"<<hx_range<<" "<<hx_bearing<<" "<<hx_range_rate <<endl;
  cout <<" msr:"<<z[0] << " "<< z[1]<< " "<<z[2]<<endl;
  VectorXd y = z - hx;

  //correct phita
  y[1] = mea_bearing - hx_bearing;
  cout <<" nor:" << mea_bearing << " y[1]" << y[1] << endl;
  // assuming there is no big jump of bearing
  if (y[1] > pi())
  {
    y[1] -= 2*pi();
  }else if(y[1] < -pi())
  {
    y[1] += 2*pi();
  }


  if(! tools.CalculateJacobian(x_,Hj_))
    return;

  MatrixXd S = Hj_ * P_ * Hj_.transpose() + Rr_;
  MatrixXd K =  P_ * Hj_.transpose() * S.inverse();

  //new estimate
  x_ = x_ + (K * y);
  unsigned int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj_) * P_;

}
