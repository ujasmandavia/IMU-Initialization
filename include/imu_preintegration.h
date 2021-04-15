#pragma onc

#include <iostream>
#include <vector>

#include <Eigen/Core>

#include "eigen_defs.h"

namespace IMU{
EIGEN_MAKE_ALIGNED_OPERATOR_NEW

const double GRAVITY_MAGNITUDE = 9.81;
const Eigen::Vector3d GRAVITY_VECTOR(0.0,0.0,-GRAVITY_MAGNITUDE);

static Eigen::Matrix6d Sigma = Eigen::Matrix6d::Identity(); // discrete

class Point{
public:
  Point() {}
  Point(const double a_x, const double a_y, const double a_z,
        const double w_x, const double w_y, const double w_z,const double dt) : a(a_x,a_y,a_z), w(w_x,w_y,w_z), dt(dt) {}
  Point(const Eigen::Vector3d& a, const Eigen::Vector3d w, const double dt) : a(a), w(w), dt(dt) {}

private:
  Eigen::Vector3d a;
  Eigen::Vector3d w;
  double dt;

}; //class Point

//Preintegration of the IMU measurements
class Preintegration{
public:
  Preintegration() {}
  Preintegration(const Eigen::Vector3d& ba, const Eigen::Vector3d& bg);
  ~Preintegration() {}
  void IntegrateNewMeasurement(const Eigen::Vector3d& a, const Eigen::Vector3d& w, const double dt);
  void SetNewGyroBias(const Eigen::Vector3d& bg);
  void SetNewAccBias(const Eigen::Vector3d& ba);
  Eigen::Vector3d GetGyroDeltaBias() const;
  Eigen::Vector3d GetGyroDeltaBias(const Eigen::Vector3d& bg) const;
  Eigen::Vector3d GetGyroOriginalBias() const;
  Eigen::Vector3d GetGyroUpdatedBias() const;
  Eigen::Vector3d GetAccDeltaBias() const;
  Eigen::Vector3d GetAccDeltaBias(const Eigen::Vector3d& ba) const;
  Eigen::Vector3d GetAccOriginalBias() const;
  Eigen::Vector3d GetAccUpdatedBias() const;
  Eigen::Matrix3d GetDeltaRotation(const Eigen::Vector3d& bg) const;
  Eigen::Vector3d GetDeltaVelocity(const Eigen::Vector3d& bg, const Eigen::Vector3d& ba) const;
  Eigen::Vector3d GetDeltaPosition(const Eigen::Vector3d& bg, const Eigen::Vector3d& ba) const;
  Eigen::Matrix3d GetUpdatedDeltaRotation() const;
  Eigen::Vector3d GetUpdatedDeltaVelocity() const;
  Eigen::Vector3d GetUpdatedDeltaPosition() const;
  Eigen::Matrix3d GetOriginalDeltaRotation() const;
  Eigen::Vector3d GetOriginalDeltaVelocity() const;
  Eigen::Vector3d GetOriginalDeltaPosition() const;

private:
  Preintegration() {}
  void Initialize(const Eigen::Vector3d& ba, const Eigen::Vector3d& bg);

public:
  double dT;
  Eigen::Matrix9d C;

  Eigen::Matrix6d Nga;

  //Values for the original bias (when integration was computed)
  Eigen::Vector6d b;
  Eigen::Matrix3d dR;
  Eigen::Vector3d dV,dP;
  Eigen::Matrix3d JRg, JVg, JVa, JPg, JPa;
  //Eigen::Vector3d avgA, avgW;
  
private:
  //Updated bias
  Eigen::Vector6d bu;
  // Dif between original and updated bias
  // This is used to compute the updated values of the preintegration
  Eigen::Vector6d db;

}; //class Preintegration

} // namespace IMU
