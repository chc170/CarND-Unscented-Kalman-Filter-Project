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
        previous_timestamp_ = meas_pkg.timestamp_;

        switch (meas_pkg.sensor_type) {
            case MeasurementPackage::RADAR:
                double rho   = meas_pkg.raw_measurements_[0];
                double phi   = meas_pkg.raw_measurements_[1];
                double rho_d = meas_pkg.raw_measurements_[2];

                VectorXd polar(3);
                polar << rho, phi, rho_d;

                VectorXd x(4);
                x << tools.PolarToCartesian(polar);

                x_ << x[0], x[1], rho_d, 0, 0;
                break;
            case MeasurementPackage::LASER:
                double px = meas_pkg.raw_measurements_[0];
                double py = meas_pkg.raw_measurements_[1];

                x_ << px, py, 0, 0, 0;
                break;
            default:
                // error
        }

        is_initialized_ = true;
        return;
    }

    // elapsed time (dt)
    float dt = (meas_pkg.timestamp_ - previous_timestamp_) * 1e-6;
    previous_timestamp_ = meas_pkg.timestamp_;

    /**
     * Prediction
     */
    if (dt > 0.001) {
        Prediction(dt);
    }

    /**
     * Update
     */
    switch (meas_pkg.sensor_type) {
        case MeasurementPackage::RADAR:
            UpdaateRadar(meas_pkg);
            break;

        case MeasurementPackage::LASER:
            UpdateLidar(meas_pkg);
            break;

        default:
            // error
    }
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
}
