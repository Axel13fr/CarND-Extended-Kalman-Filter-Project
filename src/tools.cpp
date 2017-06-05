#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

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
        assert(false);
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
    auto ro = radar_mes(RO);
    auto theta = radar_mes(THETA);
    auto state = VectorXd(STATE_SIZE);
    // x,y, Vx, Vy
    state << ro*std::cos(theta) , ro*std::sin(theta) ,
            0 , 0;

    return state;
}

Eigen::VectorXd Tools::TransformToRadarFromState(const Eigen::VectorXd &state)
{
    // state = x,y, Vx, Vy
    auto ro = sqrt(state(X)*state(X) + state(Y)*state(Y));
    auto theta = atan2(state(Y),state(X));
    assert(theta <= M_PI and theta >= -M_PI);

    auto c1 = state(X)*state(Vx) + state(Y)*state(Vy);
    assert(fabs(c1) > 0.0001);
    assert(fabs(ro) > 0.0001);
    auto ro_dot = c1 / ro;

    auto rad_mes = VectorXd(RAD_MES_SIZE);
    rad_mes << ro , theta, ro_dot;

    return rad_mes;
}

void Tools::NormalizeDeltaPhi(VectorXd& y)
{
    auto delta_phi = y(THETA);
    auto norm_delta_phi = fmod(delta_phi,2*M_PI);
    if(norm_delta_phi > M_PI){
        norm_delta_phi -= 2*M_PI;
    }else if(norm_delta_phi < -M_PI){
        norm_delta_phi += 2*M_PI;
    }else{
        // Already normalized, we are cool
    }

    assert((norm_delta_phi <= M_PI) and (norm_delta_phi >= -M_PI));
    // save norm. result !
    y(THETA) = norm_delta_phi;
}
