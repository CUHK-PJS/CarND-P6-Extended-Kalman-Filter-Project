#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter()
{
}

KalmanFilter::~KalmanFilter()
{
}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in, MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
{
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}

void KalmanFilter::Predict(double delta_T)
{
    /**
    TODO:
      * predict the state
    */

    F_(0, 0) = 1.0;
    F_(0, 1) = 0.0;
    F_(0, 2) = delta_T;
    F_(0, 3) = 0.0;

    F_(1, 0) = 0.0;
    F_(1, 1) = 1.0;
    F_(1, 2) = 0.0;
    F_(1, 3) = delta_T;

    F_(2, 0) = 0.0;
    F_(2, 1) = 0.0;
    F_(2, 2) = 1.0;
    F_(2, 3) = 0.0;

    F_(3, 0) = 0.0;
    F_(3, 1) = 0.0;
    F_(3, 2) = 0.0;
    F_(3, 3) = 1.0;

    double t1 = delta_T;
    double t2 = t1*t1;
    double t3 = t1*t2;
    double t4 = t2*t2;

    double ax = 9;
    double ay = 9;

    Q_(0, 0) = t4*ax/4.0;
    Q_(0, 1) = 0.0;
    Q_(0, 2) = t3*ax/2.0;
    Q_(0, 3) = 0.0;

    Q_(1, 0) = 0.0;
    Q_(1, 1) = t4*ay/4.0;
    Q_(1, 2) = 0.0;
    Q_(1, 3) = t3*ay/2.0;

    Q_(2, 0) = t3*ax/2.0;
    Q_(2, 1) = 0.0;
    Q_(2, 2) = t2*ax;
    Q_(2, 3) = 0.0;

    Q_(3, 0) = 0.0;
    Q_(3, 1) = t3*ay/2.0;
    Q_(3, 2) = 0.0;
    Q_(3, 3) = t2*ay;

    x_ = F_*x_;
    P_ = F_*P_*F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z)
{
    /**
    TODO:
      * update the state by using Kalman Filter equations
    */
    VectorXd y = z - H_*x_;
    MatrixXd S = H_*P_*H_.transpose() + R_;
    MatrixXd K = P_*H_.transpose()*S.inverse();
    MatrixXd I = MatrixXd::Identity(4, 4);
    x_ = x_ + K*y;
    P_ = (I-K*H_)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
    /**
    TODO:
      * update the state by using Extended Kalman Filter equations
    */
    double px = x_(0);
    double py = x_(1);
    double vx = x_(2);
    double vy = x_(3);

    VectorXd hx(3);
    hx(0) = sqrt(px*px+py*py);

    //if zero, donot change x_ and P_; to avoid nan;
    if(hx(0)<1e-6)
    {
        return;
    }

    hx(1) = atan2(py, px);
    hx(2) = (px*vx+py*vy)/hx(0);

    VectorXd y = z - hx;

    Round2PI(y(1));

    MatrixXd S = H_*P_*H_.transpose() + R_;
    MatrixXd K = P_*H_.transpose()*S.inverse();
    MatrixXd I = MatrixXd::Identity(4, 4);
    x_ = x_ + K*y;
    P_ = (I-K*H_)*P_;
}
