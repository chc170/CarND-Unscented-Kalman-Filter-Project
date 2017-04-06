#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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

  n_x_ = 5;
  n_aug = 7;
  lambda = 3 - n_aug;

  double n_sig_pts = 2 * n_aug + 1;
  weights_ = Vectord(n_sig_pts);
  Xsig_pred_ = VectorXd(n_x, n_sig_pts);

  NIS_radar_ = NULL;
  NIS_laser_ = NULL;

  previous_timestamp_ = 0;

  is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_pkg The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_pkg) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

    if (!is_initialized_) {
        Initialize(meas_pkg);
        return;
    }

    // elapsed time (dt)
    double dt = (meas_pkg.timestamp_ - previous_timestamp_) * 1e-6;
    previous_timestamp_ = meas_pkg.timestamp_;

    /**
     * Prediction
     * skip when dt is too small.
     */
    if (dt > 0.001) {
        Prediction(dt);
    }

    /**
     * Update
     */
    if (use_radar_ &&
        meas_pkg.sensor_type == MeasurementPackage::RADAR) {

        UpdaateRadar(meas_pkg);
    }
    else
    if (use_laser &&
        meas_pkg.sensor_type == MeasurementPackage::LASER) {

        UpdateLidar(meas_pkg);
    }
}

/**
 * Initialize
 * @param meas_package
 */
void UKF::Initialize(MeasurementPackage meas_pkg) {

    previous_timestamp_ = meas_pkg.timestamp_;

    auto z = meas_pkg.raw_measurements_;

    if (use_radar_ &&
        meas_pkg.sensor_type == MeasurementPackage::RADAR) {

        auto x = tools.PolarToCartesian(z);
        // px, py, rho_dot
        x_ << x(0), x(1), z(2), 0, 0;
    }
    else
    if (use_laser &&
        meas_pkg.sensor_type == MeasurementPackage::LASER) {

        // px, py
        x_ << z(0), z(1), 0, 0, 0;
    }
    else {
        return;
    }
    is_initialized_ = true;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    /**
    TODO:

    Complete this function! Estimate the object's location. Modify the state
    vector, x_. Predict sigma points, the state, and the state covariance matrix.
    */

    // Generate sigma points
    Xsig = GenerateSigmaPoints();

    // Predict sigma points
    PredictSigmaPoints(Xsig, delta_t);

    // Predict mean and covariance
    PredictMeanAndCovariance();
}

MatrixXd UKF::GenerateSigmaPoints() {

    // Augmented state matrix
    MatrixXd x_aug = VectorXd::Zero(n_aug);
    x_aug.head(n_x) = x_;

    // Augmented covariance matrix
    MatrixXd P_aug = MatrixXd::Zero(n_aug, n_aug);
    P_aug.topLeftCorner(n_x, n_x) = P_;
    P_aug(5, 5) = std_a_*std_a_;
    P_aug(6, 6) = std_yawdd_*std_yawdd_;

    // Square root matrix
    MatrixXd A_aug = MatrixXd::Zero(n_aug, n_aug);
    A_aug = P_aug.llt().matrixL();

    // Sigma points
    MatrixXd Xsig_aug = MatrixXd::Zero(n_x, 2 * n_aug + 1);
    float scale = sqrt(lambda + n_aug);
    Xsig_aug.col(0) = x_aug;

    for (int i = 0; i < n_aug; ++i) {
        Xsig_aug.col(i+1)       = x_aug + scale * A.col(i);
        Xsig_aug.col(i+1+n_aug) = x_aug - scale * A.col(i);
    }
    return Xsig_aug;
}

void UKF::PredictSigmaPoints(MatrixXd Xsig, double dt) {

    for (int i = 0; i < 2 * n_aug + 1; ++i) {

        // Extract values
        double px       = Xsig(0, i);
        double py       = Xsig(1, i);
        double v        = Xsig(2, i);
        double yaw      = Xsig(3, i);
        double yawd     = Xsig(4, i);
        double nu_a     = Xsig(5, i);
        double nu_yawdd = Xsig(6, i);

        // Predicted values
        double px_p, py_p, v_p, yaw_p, yawd_p;

        if (fabs(yawd) > 0.001) {
            px_p = px * v/yawd * (sin(yaw + yawd*dt) - sin(yaw));
            py_p = py * v/yawd * (cos(yaw) - cos(yaw + yawd*dt));
        }
        else {
            px_p = px * v * dt * cos(yaw);
            py_p = py * v * dt * sin(yaw);
        }

        v_p    = v + nu_a*dt;
        yaw_p  = yaw * yawd*dt + 0.5*nu_yawdd*dt*dt;
        yawd_p = yawd + nu_yawdd*dt;

        Xsig_pred_(0, i) = px_p;
        Xsig_pred_(1, i) = py_p;
        Xsig_pred_(2, i) = v_p;
        Xsig_pred_(3, i) = yaw_p;
        Xsig_pred_(4, i) = yawd_p;
    }
}

void UKF::PredictMeanAndCovariance() {

    // weights
    VectorXd weights = VectorXd(2*n_aug+1);
    weights(0) = lambda / (lambda+n_aug);

    for (int i = 1; i < 2*n_aug+1; ++i) {
        weights(i) = 0.5/(n_aug+lambda);
    }

    // Predict state mean
    x_.fill(0.0);
    for (int i = 0; i < 2*n_aug+1; ++i) {
        x += weights(i) * Xsig_pred_.col(i);
    }

    // Predict state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < 2*n_aug+1; ++i) {

        VectorXd diff = Xsig_pred_.col(i) - x;
        // angle normalization
        while (diff(3) >  M_PI) diff(3) -= 2. * M_PI;
        while (diff(3) < -M_PI) diff(3) += 2. * M_PI;

        P_ += weights(i) * diff * diff.transpose();
    }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_pkg
 */
void UKF::UpdateLidar(MeasurementPackage meas_pkg) {
    /**
    TODO:

    Complete this function! Use lidar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the lidar NIS.
    */

    // Transform to measurement space

    // Update
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_pkg
 */
void UKF::UpdateRadar(MeasurementPackage meas_pkg) {
    /**
    TODO:

    Complete this function! Use radar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the radar NIS.
    */

    // Transform to measurement space

    // Update
}
