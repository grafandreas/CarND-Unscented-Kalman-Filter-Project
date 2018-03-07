#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

#define FLOAT_FIX(x) ((fabs(x)) < 0.00001 ? 0.00001 : (x))
#define NORMALIZE_PI(x) while ((x) > M_PI) {\
  x -= 2. * M_PI;\
}\
  while ((x) < -M_PI) {\
  x += 2. * M_PI;\
}\

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

  VectorXd polar2Cart(const VectorXd& polarC);

};

#endif /* TOOLS_H_ */
