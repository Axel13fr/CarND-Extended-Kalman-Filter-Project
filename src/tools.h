#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
    enum state_E{
        X = 0,
        Y,
        Vx,
        Vy,

        STATE_SIZE
    };
    enum rad_mes_E{
      RO = 0,
      THETA,
      RO_DOT,

      RAD_MES_SIZE
    };

  /**
  * A helper method to calculate RMSE.
  */
  static VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  static MatrixXd CalculateJacobian(const VectorXd& x_state);

  /**
   * @brief CalculatePosFromRadar
   * @param radar_mes see enum for content
   * @return state vector init using radar position
   */
  static VectorXd CalculatePosFromRadar(const VectorXd& radar_mes);

  /**
   * @brief Tools::CalculateRadarMesFromState
   * @param state see enum for content
   * @return radar measurement based on state vector
   */
  static Eigen::VectorXd TransformToRadarFromState(const Eigen::VectorXd &state);

  static void NormalizeDeltaPhi(Eigen::VectorXd &y);
};

#endif /* TOOLS_H_ */
