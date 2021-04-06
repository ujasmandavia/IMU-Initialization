#include "../include/imu_preintegration.h"
#include "../include/so3.h"

namespace IMU{

//Integration of 1 gyro measurement
class IntegratedRotation{
public:
  IntegratedRotation() = delete;
  IntegratedRotation(const Eigen::Vector3d& w, const Eigen::Vector3d& bias, const double time) : deltaT(time){
    //The below lines are the preintegrated rotation measurement
    const Eigen::Vector3d dr = time*(w-bias); 
    deltaR = ExpSO3(dr.x(),dr.y(),dr.z());  //Eq 28 in the main preintegration paper
    rightJ = RightJacobianSO3(dr.x(),dr.y(),dr.z());
  }

public:
  const double deltaT;  //integration time
  Eigen::Matrix3d deltaR;  //integrated rotation
  Eigen::Matrix3d rightJ;  //right jacobian
}; //class IntegratedRotation

Preintegration::Preintegration(const Eigen::Vector3d& ba, const Eigen::Vector3d& bg){
  Initialize(bg,bg);
}

void Preintegration::IntegrateNewMeasurement(const Eigen::Vector3d& a, const Eigen::Vector3d& w, const double dt){
  //Position is updated first since it depends on previously computer velocity and rotation
  //velocity is updated second, as it depends on previously computed rotation
  //Rotation is the last to be updated

  //Matrices to compute covariances
  Eigen::Matrix9d A;
  A.setIdentity();

}

} //namespace IMU
