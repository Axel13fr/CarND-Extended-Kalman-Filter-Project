#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

    auto rmse = VectorXd(4);
    rmse << 0,0,0,0;
    if(estimations.size() != ground_truth.size() or estimations.size() == 0){
        std::cout << "size error in RMSE calc " << std::endl;
        return rmse;
    }

    for(unsigned int i=0; i < estimations.size() ; i++){
        VectorXd residual = estimations[i] - ground_truth[i];
        residual = residual.array()*residual.array();
        rmse += residual;
    }

    rmse = rmse/estimations.size();
    rmse = rmse.array().sqrt();

    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    MatrixXd Hj(3,4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    //pre-compute a set of terms to avoid repeated calculation
    float c1 = px*px+py*py;
    float c2 = sqrt(c1);
    float c3 = (c1*c2);

    //check division by zero
    if(fabs(c1) < 0.0001){
        std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
        return Hj;
    }

    //compute the Jacobian matrix
    Hj << (px/c2), (py/c2), 0, 0,
          -(py/c1), (px/c1), 0, 0,
          py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

    return Hj;
}

Eigen::VectorXd Tools::CalculatePosFromRadar(const Eigen::VectorXd &radar_mes)
{
    auto ro = radar_mes(0);
    auto theta = radar_mes(1);
    auto ro_dot = radar_mes(2);
    auto state = VectorXd(4);
    // x,y, Vx, Vy
    state << ro*std::cos(theta) , ro*std::sin(theta) ,
            0 , 0;

    return state;
}


