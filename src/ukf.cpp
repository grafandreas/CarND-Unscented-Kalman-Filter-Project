#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {

  previous_timestamp_ = 0.0;
  time_us_ = 0;

  is_initialized_ = false;

  // Parameterization of behavior, use or not use laser / radar
  use_laser_ = true;
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;  // A deviation of 30 for acc would be a huge problem, it would be way to imprecise

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.60;  // Similar, this would be a much to high value

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  // Calculation of important variables just as in project lesson
  n_x_ = x_.size();
  n_aug_ = n_x_ + 2;
  n_sig_points_ = 2 * n_aug_ + 1;
  lambda_ = 3 - n_aug_;

  // Note that setting zero is very important
  Xsig_pred_ = MatrixXd(n_x_, n_sig_points_);
  Xsig_pred_.setZero();


  // Weights of sigma points
  weights_ = VectorXd(n_sig_points_);
  weights_.setZero();

  // covariance matrices initialization
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_ * std_radr_, 0, 0, 0, std_radphi_ * std_radphi_, 0, 0, 0, std_radrd_
      * std_radrd_;
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_ * std_laspx_, 0, 0, std_laspy_ * std_laspy_;

}

UKF::~UKF() {
}

// This initializes all our parameters. Called only once on startup.

void UKF::Init(MeasurementPackage measurement_pack) {
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    Tools t;

    // polar2Cart calculates All, but we only use X/Y
    //
    VectorXd c = t.polar2Cart(measurement_pack.raw_measurements_);
    x_ = c;
  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    x_ << measurement_pack.raw_measurements_[0], measurement_pack
        .raw_measurements_[1], 0.0, 0.0,  // Should 2 and 3 should be 0 anyway, but use constant so
    0.0;
  }

  // The first time_stamp is being set
  //
  previous_timestamp_ = measurement_pack.timestamp_;

  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < weights_.size() ; i++)
    weights_(i) = 0.5 / (n_aug_ + lambda_);


  // done initializing, no need to predict or update
  is_initialized_ = true;
//  cout << "Init" <<endl;
  return;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
// cout << " -> PM" <<endl;

  if (!is_initialized_) {
    Init(measurement_pack);
  }

  auto delta_t = ((measurement_pack.timestamp_ - previous_timestamp_)
      / 1000000.0);
  previous_timestamp_ = measurement_pack.timestamp_;

  Prediction(delta_t);

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(measurement_pack);
  } else if(measurement_pack.sensor_type_ == MeasurementPackage::LASER &&  use_laser_) {
    UpdateLidar(measurement_pack);
  }
  assert(!x_.hasNaN());

}

// This updates the covariance matrix.
//
void UKF::UpdateCovarianceMatrix() {
  P_.setZero();
  for (int i = 0; i < n_sig_points_; i++) {  //iterate over sigma points
    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // Angle normalization
    NormAng(&(x_diff(3)));
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}


void UKF::PredictSigmaPoints(double delta_t) {
  double delta_t2 = delta_t * delta_t;

  VectorXd x_aug = VectorXd(n_aug_);

   MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

   MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_points_);


   x_aug.fill(0.0);
   x_aug.head(n_x_) = x_;


   P_aug.fill(0);
   P_aug.topLeftCorner(n_x_, n_x_) = P_;
   P_aug(5, 5) = std_a_ * std_a_;
   P_aug(6, 6) = std_yawdd_ * std_yawdd_;


   // Square root of P matrix
   MatrixXd L = P_aug.llt().matrixL();
   // Create sigma points
   Xsig_aug.col(0) = x_aug;
   double sqrt_lambda_n_aug = sqrt(lambda_ + n_aug_);  // Save some computations
   VectorXd sqrt_lambda_n_aug_L;
   for (int i = 0; i < n_aug_; i++) {
     sqrt_lambda_n_aug_L = sqrt_lambda_n_aug * L.col(i);
     Xsig_aug.col(i + 1) = x_aug + sqrt_lambda_n_aug_L;
     Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt_lambda_n_aug_L;
   }


  for (int i = 0; i < n_sig_points_; i++) {

    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);
    // Precalculate sin and cos
    double sin_yaw = sin(yaw);
    double cos_yaw = cos(yaw);
    double arg = yaw + yawd * delta_t;


    double px_p, py_p;
    //numeric trick for div by zero
    if (fabs(yawd) > 0.00001) {
      double v_yawd = v / yawd;
      px_p = p_x + v_yawd * (sin(arg) - sin_yaw);
      py_p = p_y + v_yawd * (cos_yaw - cos(arg));
    } else {
      double v_delta_t = v * delta_t;
      px_p = p_x + v_delta_t * cos_yaw;
      py_p = p_y + v_delta_t * sin_yaw;
    }
    double v_p = v;
    double yaw_p = arg;
    double yawd_p = yawd;

    // noise
    px_p += 0.5 * nu_a * delta_t2 * cos_yaw;
    py_p += 0.5 * nu_a * delta_t2 * sin_yaw;
    v_p += nu_a * delta_t;
    yaw_p += 0.5 * nu_yawdd * delta_t2;
    yawd_p += nu_yawdd * delta_t;

    // Write predicted sigma point
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }
}

void UKF::Prediction(double delta_t) {
//  cout << "-> PR" <<endl;



  PredictSigmaPoints(delta_t);



  // Predicted state mean
  x_ = Xsig_pred_ * weights_;
  assert(!x_.hasNaN());
  // Predicted state covariance matrix
  UpdateCovarianceMatrix();


//  cout << "<-PR" <<endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
//  cout << "UR" <<endl;
  // measurement dimension, radar : r, phi, and r_dot
  int n_z = 3;
  //  sigma point matrix in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sig_points_);
  // Transform sig points to measurement space
  assert(!Xsig_pred_.hasNaN());
  for (int i = 0; i < n_sig_points_; i++) {

    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);

    p_x = FLOAT_FIX(p_x);
    p_y = FLOAT_FIX(p_y);

    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);          //r
    Zsig(1, i) = atan2(p_y, p_x);                   //phi
    Zsig(2, i) = (p_x * v1 + p_y * v2) / Zsig(0, i);   //r_dot

    // Fail early in case of calculation problems
    //
    if(isnan(Zsig(2,i))) {
      cout << "ZsigNan " << p_x << ","<< p_y << "," << v1 << "," << v2 << "," << Zsig(0, i) << endl;
    }

  }
//  cout << "zsig" << Zsig << endl;
  UpdateUKF(meas_package, Zsig, n_z);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
//  cout << "UL" <<endl;
  // measurement dimension
  int n_z = 2;

  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z, n_sig_points_);
  UpdateUKF(meas_package, Zsig, n_z);
}

// Universal update function
void UKF::UpdateUKF(MeasurementPackage measure_arg, MatrixXd Zsig, int measurement_dim) {
//  cout << "UU" <<endl;
  assert(!Zsig.hasNaN());

  // Mean predicted measurement
  VectorXd z_pred = VectorXd(measurement_dim);
  z_pred = Zsig * weights_;
  //measurement covariance matrix
  MatrixXd S = MatrixXd(measurement_dim, measurement_dim);
  S.fill(0.0);
  for (int i = 0; i < n_sig_points_; i++) {

    VectorXd z_diff = Zsig.col(i) - z_pred;
    // Normalize angle
    NormAng(&(z_diff(1)));
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  // Add measurement noise covariance matrix
  MatrixXd R = MatrixXd(measurement_dim, measurement_dim);
  if (measure_arg.sensor_type_ == MeasurementPackage::RADAR) {  // Radar
//    cout << "Radar" << endl;
    R = R_radar_;
  } else if (measure_arg.sensor_type_ == MeasurementPackage::LASER) {  // Lidar
//    cout << "Laser" << endl;
    R = R_laser_;
  }
//  cout << S.cols() << " " << S.rows() << endl;
//  cout << R.cols() << " " << R.rows() << endl;

  S = S + R;

  // Create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, measurement_dim);
  //  cross correlation matrix calculation
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_points_; i++) {

    VectorXd z_diff = Zsig.col(i) - z_pred;
    if (measure_arg.sensor_type_ == MeasurementPackage::RADAR) {  // Radar
      // Angle normalization
      NormAng(&(z_diff(1)));
    }
    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // Angle normalization
    NormAng(&(x_diff(3)));
    assert(!weights_.hasNaN());
    assert(!x_diff.hasNaN());
    assert(!z_diff.hasNaN());
    assert(!z_diff.transpose().hasNaN());

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  VectorXd z = measure_arg.raw_measurements_;
  //gain K;
  assert(!Tc.hasNaN());
  assert(!S.inverse().hasNaN())
;  MatrixXd K = Tc * S.inverse();
  //
  VectorXd z_diff = z - z_pred;
  if (measure_arg.sensor_type_ == MeasurementPackage::RADAR) {  // Radar
    // Angle normalization
    NormAng(&(z_diff(1)));
  }
  //  state mean and covariance matrix update
  assert(!K.hasNaN());
  assert(!z_diff.hasNaN());
  x_ = x_ + K * z_diff;
  assert(!x_.hasNaN());
  P_ = P_ - K * S * K.transpose();
  // NIS calculation
  if (measure_arg.sensor_type_ == MeasurementPackage::RADAR) {  // Radar
    NIS_radar_ = z.transpose() * S.inverse() * z;
  } else if (measure_arg.sensor_type_ == MeasurementPackage::LASER) {  // Lidar
    NIS_laser_ = z.transpose() * S.inverse() * z;
  }
}

// Normalize an Angle. This could be changed to an inline function for
// production.
//
void UKF::NormAng(double *ang) {

  // Contains some sanity checks because of some original problems
  // that came from incorrect initalization.
  //
  if(*ang > 500) {
    cout << "Angle Normalization problem : " << *ang << endl;
    *ang = 500;
  }

  if(*ang < -500) {
    cout << "Angle Normalization problem : " << *ang << endl;
    *ang = -500;
  }
  while (*ang > M_PI) {
//    cout << * ang << endl;
    *ang -= 2. * M_PI;
  }
  while (*ang < -M_PI)
    *ang += 2. * M_PI;
}


