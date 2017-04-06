#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size

    //accumulate squared residuals
    for(int i=0; i < estimations.size(); ++i){
        VectorXd residual = ground_truth[i]-estimations[i];
        residual = residual.array() * residual.array();

        rmse += residual;
    }

    //calculate the mean
    rmse = rmse.array() / estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();
    return rmse;
}

MatrixXd Tools::PolarToCartesian(const VectorXd& x) {

    // recover state parameters
    float rho = x[0];
    float phi = x[1];
    float rho_dot = x[2];

    // coordinate conversion calculations
    float px = rho * cos(phi);
    float py = rho * sin(phi);

    float vx = rho_dot * cos(phi);
    float vy = rho_dot * sin(phi);

    VectorXd x_c(4);
    x_c << px, py, vx, vy;

    return x_c;
}
