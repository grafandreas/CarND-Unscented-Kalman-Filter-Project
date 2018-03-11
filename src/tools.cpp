#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
 // Code taken from pervious project

//  cout << "-> RMSE " << estimations.size() << endl;



  assert(estimations.size() == ground_truth.size());

  VectorXd rmse = VectorXd(4);
  rmse << 0.0,0.0,0.0,0.0;

//  cout << "I " << estimations.size() << endl;
//  cout << estimations.back() << endl;
//  cout << "----" << endl;
//  cout  << ground_truth.back() << endl;
//  cout << "++++" << endl;

  for(int i = 0; i < estimations.size(); i++) {


    VectorXd delta = estimations[i] - ground_truth[i];
    VectorXd sq = delta.array() * delta.array();

    rmse+=sq;
  }
//  cin.get();
  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();

//  assert(rmse(0) < 0.5);
//  assert(rmse(1) < 0.5);
//  cout << "<- RMSE" << rmse << endl;
  return rmse;

}


// Taken from the previous project
VectorXd Tools::polar2Cart(const VectorXd& polarC) {
  VectorXd res = VectorXd(5);


  auto vx =     polarC(2) * cos(polarC(1));
  auto vy =     polarC(2) * sin(polarC(1));



  auto v =  sqrt(vx*vx+vy*vy);


  res << polarC(0) * cos(polarC(1)), // x
      polarC(0) * sin(polarC(1)),    // y
      v,
      0.0,
      0.0;
  return res;
}


