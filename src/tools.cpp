#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools()
{
}

Tools::~Tools()
{
}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth)
{
    /**
    TODO:
      * Calculate the RMSE here.
    */
    VectorXd rmse(4);
    rmse << 0.0, 0.0, 0.0, 0.0;

    if (estimations.size() == 0 || estimations.size() != ground_truth.size())
    {
        std::cout << "invalid estimation or ground truth data" << std::endl;
        return rmse;
    }

    for (int i = 0; i < estimations.size(); ++i)
    {
        VectorXd residual = estimations[i] - ground_truth[i];
        rmse.array() += residual.array() * residual.array();;
    }
    rmse /= estimations.size();
    rmse = rmse.array().sqrt();

    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state)
{
    /**
    TODO:
      * Calculate a Jacobian here.
    */
    MatrixXd Hj(3, 4);

    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);

    double c1 = std::max(px * px + py * py, 1e-6);//to avoid nan;
    double c2 = sqrt(c1);
    double c3 = c1*c2;

    Hj(0, 0) = px/c2;
    Hj(0, 1) = py/c2;
    Hj(0, 2) = 0.0;
    Hj(0, 3) = 0.0;

    Hj(1, 0) = -py/c1;
    Hj(1, 1) = px/c1;
    Hj(1, 2) = 0.0;
    Hj(1, 3) = 0.0;

    Hj(2, 0) = py*(vx*py-vy*px)/c3;
    Hj(2, 1) = px*(vy*px-vx*py)/c3;
    Hj(2, 2) = px/c2;
    Hj(2, 3) = py/c2;

    return Hj;
}
