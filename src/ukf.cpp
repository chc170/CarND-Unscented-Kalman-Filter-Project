#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

static const double twoPI = 2 * M_PI;

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
    P_.fill(0.0);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 3.;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = M_PI / 2;

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

    n_x_    = 5;
    n_aug_  = 7;
    n_sig_  = 2 * n_aug_ + 1;
    lambda_ = 3 - n_aug_;

    weights_   = VectorXd::Zero(n_sig_);
    Xsig_pred_ = MatrixXd::Zero(n_x_, n_sig_);

    NIS_radar_ = 0;
    NIS_laser_ = 0;

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
     */

    //while (dt > 0.1) {
    //    Prediction(0.05);
    //    dt -= 0.05;
    //}
    Prediction(dt);

    /**
     * Update
     */
    if (use_radar_ &&
        meas_pkg.sensor_type_ == MeasurementPackage::RADAR) {

        UpdateRadar(meas_pkg);
    }
    else
    if (use_laser_ &&
        meas_pkg.sensor_type_ == MeasurementPackage::LASER) {

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
        meas_pkg.sensor_type_ == MeasurementPackage::RADAR) {

        auto x = tools.PolarToCartesian(z);

        // px, py, rho_dot
        x_ << x(0), x(1), z(2), 0, 0;
    }
    else
    if (use_laser_ &&
        meas_pkg.sensor_type_ == MeasurementPackage::LASER) {

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
    MatrixXd Xsig = GenerateSigmaPoints();

    // Predict sigma points
    PredictSigmaPoints(Xsig, delta_t);

    // Predict mean and covariance
    PredictMeanAndCovariance();
}

MatrixXd UKF::GenerateSigmaPoints() {

    // Augmented state matrix
    VectorXd x_aug = VectorXd::Zero(n_aug_);
    x_aug.head(n_x_) = x_;

    // Augmented covariance matrix
    MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(5, 5) = std_a_*std_a_;
    P_aug(6, 6) = std_yawdd_*std_yawdd_;

    // Square root matrix
    MatrixXd A_aug = P_aug.llt().matrixL();

    // Sigma points
    MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, n_sig_);
    double scale = sqrt(lambda_ + n_aug_);
    Xsig_aug.col(0) = x_aug;

    for (int i = 0; i < n_aug_; ++i) {
        Xsig_aug.col(i+1)        = x_aug + scale * A_aug.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - scale * A_aug.col(i);
    }

    return Xsig_aug;
}

void UKF::PredictSigmaPoints(MatrixXd Xsig, double dt) {

    for (int i = 0; i < n_sig_; ++i) {

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
            px_p = px + v/yawd * (sin(yaw + yawd*dt) - sin(yaw));
            py_p = py + v/yawd * (cos(yaw) - cos(yaw + yawd*dt));
        }
        else {
            px_p = px + v * dt * cos(yaw);
            py_p = py + v * dt * sin(yaw);
        }

        px_p += 0.5*nu_a*dt*dt * cos(yaw);
        py_p += 0.5*nu_a*dt*dt * sin(yaw);

        v_p    = v + nu_a*dt;
        yaw_p  = yaw + yawd*dt + 0.5*nu_yawdd*dt*dt;
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
    weights_(0) = lambda_ / (lambda_+n_aug_);

    for (int i = 1; i < n_sig_; ++i) {
        weights_(i) = 0.5/(n_aug_+lambda_);
    }

    // Predict state mean
    x_.fill(0.0);
    for (int i = 0; i < n_sig_; ++i) {
        x_ += weights_(i) * Xsig_pred_.col(i);
    }

    // Predict state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < n_sig_; ++i) {
        VectorXd x_diff = Xsig_pred_.col(i) - x_;

        // angle normalization
        x_diff(3) += M_PI;
        x_diff(3) = x_diff(3) - twoPI * floor(x_diff(3)/twoPI) - M_PI;

        P_ += weights_(i) * x_diff * x_diff.transpose();
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
    int n_z = 2;

    // Transform sigma points to measurement space
    MatrixXd Zsig = MatrixXd(n_z, n_sig_);
    for (int i = 0; i < n_sig_; ++i) {

        double px  = Xsig_pred_(0, i);
        double py  = Xsig_pred_(1, i);

        // px, py
        Zsig(0, i) = px;
        Zsig(1, i) = py;
    }

    // mean predicted measurement
    VectorXd z_pred = VectorXd::Zero(n_z);
    for (int i = 0; i < n_sig_; ++i) {
        z_pred += weights_(i) * Zsig.col(i);
    }

    // measurement covariance matrix S
    MatrixXd S = MatrixXd::Zero(n_z, n_z);
    for (int i = 0; i < n_sig_; ++i) {
        // residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        // angle normalization
        z_diff(1) += M_PI;
        z_diff(1) = z_diff(1) - twoPI * floor(z_diff(1)/twoPI) - M_PI;

        S += weights_(i) * z_diff * z_diff.transpose();
    }

    // measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z, n_z);
    R << std_laspx_*std_laspx_, 0,
         0, std_laspy_*std_laspy_;
    S += R;

    // Update
    // correlation matrix
    MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
    for (int i = 0; i < n_sig_; ++i) {

        // residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        // angle normalization
        z_diff(1) += M_PI;
        z_diff(1) = z_diff(1) - twoPI * floor(z_diff(1)/twoPI) - M_PI;

        // state diff
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        // angle normalization
        x_diff(3) += M_PI;
        x_diff(3) = x_diff(3) - twoPI * floor(x_diff(3)/twoPI) - M_PI;

        Tc += weights_(i) * x_diff * z_diff.transpose();
    }

    // Kalman gain K
    MatrixXd K = Tc * S.inverse();

    // residual
    auto z = meas_pkg.raw_measurements_;
    VectorXd z_diff = z - z_pred;

    // angle normalization
    z_diff(1) += M_PI;
    z_diff(1) = z_diff(1) - twoPI * floor(z_diff(1)/twoPI) - M_PI;

    // update state mean and covariance matrix
    x_ += K * z_diff;
    P_ -= K * S * K.transpose();

    // NIS
    NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
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
    int n_z = 3;

    // Transform sigma points to measurement space
    MatrixXd Zsig = MatrixXd(n_z, n_sig_);
    for (int i = 0; i < n_sig_; ++i) {

        double px  = Xsig_pred_(0, i);
        double py  = Xsig_pred_(1, i);
        double v   = Xsig_pred_(2, i);
        double yaw = Xsig_pred_(3, i);

        double v1 = cos(yaw) * v;
        double v2 = sin(yaw) * v;

        // rho, phi, rho_dot
        double rho = sqrt(px*px + py*py);
        double rho_dot = 0;
        if (rho != 0) {
            rho_dot = (px*v1 + py*v2) / rho;
        }
        Zsig(0, i) = rho;
        Zsig(1, i) = atan2(py, px);
        Zsig(2, i) = rho_dot;
    }

    // mean predicted measurement
    VectorXd z_pred = VectorXd::Zero(n_z);
    for (int i = 0; i < n_sig_; ++i) {
        z_pred += weights_(i) * Zsig.col(i);
    }

    // measurement covariance matrix S
    MatrixXd S = MatrixXd::Zero(n_z, n_z);
    for (int i = 0; i < n_sig_; ++i) {
        // residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        // angle normalization
        z_diff(1) += M_PI;
        z_diff(1) = z_diff(1) - twoPI * floor(z_diff(1)/twoPI) - M_PI;

        S += weights_(i) * z_diff * z_diff.transpose();
    }

    // measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z, n_z);
    R << std_radr_*std_radr_, 0, 0,
         0, std_radphi_*std_radphi_, 0,
         0, 0, std_radrd_*std_radrd_;
    S += R;

    // Update
    // correlation matrix
    MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
    for (int i = 0; i < n_sig_; ++i) {

        // residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        // angle normalization
        z_diff(1) += M_PI;
        z_diff(1) = z_diff(1) - twoPI * floor(z_diff(1)/twoPI) - M_PI;

        // state diff
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        // angle normalization
        x_diff(3) += M_PI;
        x_diff(3) = x_diff(3) - twoPI * floor(x_diff(3)/twoPI) - M_PI;

        Tc += weights_(i) * x_diff * z_diff.transpose();
    }

    // Kalman gain K
    MatrixXd K = Tc * S.inverse();

    // residual
    auto z = meas_pkg.raw_measurements_;
    VectorXd z_diff = z - z_pred;

    // angle normalization
    z_diff(1) += M_PI;
    z_diff(1) = z_diff(1) - twoPI * floor(z_diff(1)/twoPI) - M_PI;

    // update state mean and covariance matrix
    x_ += K * z_diff;
    P_ -= K * S * K.transpose();

    // NIS
    NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}
