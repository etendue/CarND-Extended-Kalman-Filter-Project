#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  bool CalculateJacobian(const VectorXd& x_state,MatrixXd &Hj);
  /**
   * A help method to normalize a cyclic value e.g. angle in -pi ~ pi
   */
  double NormalizeMinMax(double x, double min, double max);

};

#endif /* TOOLS_H_ */
