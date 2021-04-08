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
  Initialize(ba,bg);
}

void Preintegration::IntegrateNewMeasurement(const Eigen::Vector3d& a, const Eigen::Vector3d& w, const double dt){
  //Position is updated first since it depends on previously computer velocity and rotation
  //velocity is updated second, as it depends on previously computed rotation
  //Rotation is the last to be updated

  //Matrices to compute covariances
  Eigen::Matrix9d A;
  A.setIdentity();

  Eigen::Matrix<double,9,6> B;
  B.setZero();

  Eigen::Vector3d acc = a - b.tail<3>();

  //Update the delta position dP and velocity dV (rely on no-updated delta rotation)
  dP += dV*dt + 0.5*dR*acc*dt*dt;
  dV += dR*acc*dt;

  //Compute velocity and position part of matrices A and B(rely on non-updated rotation)
  Eigen::Matrix3d Wacc = Skew(acc);
  A.block<3,3>(3,0) = -dR*dt*Wacc;
  A.block<3,3>(6,0) = -0.5*dR*dt*dt*Wacc;
  A.block<3,3>(6,3) = Eigen::Matrix3d::Identity()*dt;
  B.block<3,3>(3,3) = dR*dt;
  B.block<3,3>(6,3) = 0.5*dR*dt*dt;

  //Update the position and velocity jacobian wrt bias updates
  JPa = JPa + JVa*dt -0.5*dR*dt*dt;
  JPg = JPg + JVg*dt -0.5*dR*dt*dt*Wacc*JRg;
  JVa = JVa - dR*dt;
  JVg = JVg - dR*dt*Wacc*JRg;

  //Update delta rotation
  IntegratedRotation dRi(w,b.head<3>(),dt);
  dR *= dRi.deltaR;

  //Compute rotation part of martices A and B
  A.block<3,3>(0,0) = dRi.deltaR.transpose();
  B.block<3,3>(0,0) = dRi.rightJ*dt;

  //Update covariance
  C = A*C*A.transpose() + B*Sigma*B.transpose();

  //Update rotation jacobian wrt bias correction
  JRg = dRi.deltaR.transpose()*JRg - dRi.rightJ*dt;

  //Total integrated time
  dT += dt;
}

void Preintegration::SetNewGyroBias(const Eigen::Vector3d& bg){
  bu.head<3>() = bg;
  db.head<3>() = bg - b.head<3>();
}

void Preintegration::SetNewAccBias(const Eigen::Vector3d& ba){
  bu.tail<3>() = ba;
  db.tail<3>() = ba - b.tail<3>();
}

Eigen::Vector3d Preintegration::GetGyroDeltaBias() const{
  return db.head<3>();
}

Eigen::Vector3d Preintegration::GetGyroDeltaBias(const Eigen::Vector3d& bg) const{
  return bg - b.head<3>();
}

Eigen::Vector3d Preintegration::GetGyroOriginalBias() const{
  return b.head<3>();
}

Eigen::Vector3d Preintegration::GetGyroUpdatedBias() const{
  return bu.head<3>();
}

Eigen::Vector3d Preintegration::GetAccDeltaBias() const{
  return db.tail<3>()
}

Eigen::Vector3d Preintegration::GetAccDeltaBias(const Eigen::Vector3d& ba) const{
  return ba - b.tail<3>();
}

Eigen::Vector3d Preintegration::GetAccOriginalBias() const{
  return b.tail<3>();
}

Eigen::Vector3d Preintegration::GetAccUpdatedBias() const{
  return bu.tail<3>();
}

Eigen::Matrix3d Preintegration::GetDeltaRotation(const Eigen::Vector3d& bg) const{
  return dR*ExpSO3(JRg*(bg-b.head<3>()));  //equ A.14
}

Eigen::Vector3d Preintegration::GetDeltaVelocity(const Eigen::Vector3d& bg, const Eigen::Vector3d& ba) const{
  return dV + JVg*(bg-b.head<3>()) + JVa*(ba-b.tail<3>());  //equ A.18
}

Eigen::Vector3d Preintegration::GetDeltaPosition(const Eigen::Vector3d& bg, const Eigen::Vector3d& ba) const{
  return dP + JPg*(bg - b.head<3>()) + JPa*(ba - b.tail<3>());  //equ A.19
}

Eigen::Matrix3d Preintegration::GetUpdatedDeltaRotation() const{
  return dR*ExpSO3(JRg*db.head<3>());
}

Eigen::Vector3d Preintegration::GetUpdatedDeltaVelocity() const{
  return dV + JVg*db.head<3>() + JVa*db.tail<3>();
}

Eigen::Vector3d Preintegration::GetUpdatedDeltaPosition() const{
  return dP + JPg*db.head<3>() + JPa*db.tail<3>();
}

Eigen::Matrix3d Preintegration::GetOriginalDeltaRotation() const{
  return dR;
}

Eigen::Vector3d Preintegration::GetOriginalDeltaVelocity() const{
  return dV;
}

Eigen::Vector3d Preintegration::GetOriginalDeltaPosition() const{
  return dP;
}

void Preintegration::Initialize(const Eigen::Vector3d& ba, const Eigen::Vector3d& ba){
  dT = 0.0;
  C.setZero();

  b.head<3>() = bg;
  b.tail<3>() = ba;
  dR.setIdentity();
  dV.setZero();
  dP.setZero();
  JRg.setZero();
  JVg.setZero();
  JVa.setZero();
  JPg.setZero();
  JPa.setZero();
  //avgA.setZero();
  //avgW.setZero();

  bu = b;
  db.setZero();
}

} //namespace IMU
