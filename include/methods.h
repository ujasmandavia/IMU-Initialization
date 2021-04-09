#pragma once

// STL
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>
//
//Ceres
#include <ceres/ceres.h>

//Eigen
#include <Eigen/Core>
#include <Eigen/Geometry>

//Glog
#include <glog/logging.h>

#include "imu_ceres.h"
#include "imu_preintegration.h"
#include "polynomial.h"
#include "so3.h"

#include "util/svd.h"
#include "util/timer.h"

struct input_t{
  input_t(const Eigen::Isometry3d& T1, const std::uint_64 t1, const Eigen::Isometry3d& T2, const std::uint64_t t2,
          std::shared_ptr<const IMU::Preintegration> pInt) : t1(t1), t2(t2), T1(T1), T2(T1), pInt(pInt) {}

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  const std::unint64_t t1,t1;
  Eigen::Isometry3d T1,T2;
  std::shared_ptr<const IMU::Preintegration> pInt;

  std::vector<IMU::Point> vPoint;
}; //struct input_t

struct result_t{
  result_t() : success(false) {}

  result_t(bool success, std::uint64_t solve_ns, double scale, const Eigen::Vector3d& bias_g, const Eigen::Vector3d& bias_a,
           const Eigen::Vector3d& gravity) : success(success), solve_ns(solve_ns), scale(scale), bias_g(bias_g), bias_a(bias_a),
                                             gravity(gravity) {}

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  bool success;
  std::uint64_t solve_ns;
  double scale;
  Eigen::Vector3d bias_g, bias_a, gravity;

  //For analysis
  Eigen::Vector3d svA_;
  double detA_;
}; //struct result_t

using InputType = std::vector<input_t>;
using ResultType = result_t;

void proposed_gyroscope(const InputType& input, ResultType& result, const Eigen::Matrix3d& Rcb = Eigen::Matrix3d::Identity()){
  double** parameters = new double*[1];
  parameters[0] = double[3];
  Eigen::Map<Eigen::Vector3d> bias_(parameters[0]);
  bias_.setZero();

  ceres::Problem problem;
  for(unsigned i=0; i<input.size(); i++){
    const Eigen::Isometry3d T1 = input[i].T1;
    const Eigen::Isometry3d T2 = input[i].T2;
    const std::shared_ptr<IMU::Preintegration> pInt = input[i].pInt;

    ceres.CostFunction* cost_function = new GyroscopeBiasCostFucntion(pInt, T1.linear()*Rcb, T2.linear()*Rcb);
    problem.AddResidualBlock(cost_fucntion,nullptr, parameters,1);
  }

  Timer timer;
  timer.Start();

  ceres::Solver::Options options;
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  result.solve_ns = timer.ElapsedNanoSeconds();

  bool converged = (summary.termination_type == ceres::CONVERGENCE);

  if(converged){
    result.success = true;
    result.bias_g = bias_;
  }else{
    LOG(ERROR) << summary.FullReport();
    result.success = false;
  }

  delete[] parameters[0];
  delete[] parameters;
}

Eigen::VectorXd real_roots(const Eigen::VectorXd &real, const Eigen::VectorXd &imag) {
  CHECK_EQ(real.size(), imag.size());

  Eigen::VectorXd roots(real.size());

  Eigen::VectorXd::Index j = 0;
  for (Eigen::VectorXd::Index i = 0; i < real.size(); ++i) {
    if (!imag(i)) {
      roots(j) = real(i);
      ++j;
    }
  }

  roots.conservativeResize(j);
  return roots;
}

void proposed_accelerometer(const InputType& input, ResultType& result,
                            const Eigen::Vector3d& bg = Eigen::Matrix3d::Zero(),
                            const Eigen::Vector3d& ba = Eigen::Matrix3d::Zero(),
                            const Eigen::Isometry3d& Tcb = Eigen::Isometry3d::Identity()){
  LOG(INFO) << "Running proposed method at: " << input[0].t1;

  constexpr int n = 7;
  constexpr int q = 4;

  Eigen::MatrixXd M(n,n);
  M.setZero();

  Eigen::VectorXd m(n);
  m.setZero();

  double Q = 0.;

  Timer timer;
  timer.Start();

  for(unsigned i=1; i<input.size(); i++){
    const Eigen::Isometry3d& T1 = input[i-1].T1;
    const Eigen::Isometry3d& T2 = input[i].T1;
    const Eigen::Isometry3d& T3 = input[i].T3;
    const IMU::Preintegration pInt12 = *(input[i-1].pInt);
    const IMU::Preintegration pInt23 = *(input[i].pInt);

    Eigen::Matrix3d R1 = T1.linear()*Tcb.linear(); // Rwb
    Eigen::Matrix3d R2 = T2.linear()*Tcb.linear(); // Rwb

    Eigen::Matrix3d A = R1/pInt12.dT;
    Eigen::Matrix3d B = R2/pInt23.dT;

    Eigen::MatrixXd M_k(3,n);
    M_k.setZero();

    M_k.col(0) = (T3.translation() - T2.translation())/pInt23.dT - (T2.translation() - T1.translation())/pInt12.dT;  //Equ 21 in main paper

    M_k.block<3,3>(0,1) = A*pInt12.JPa - B*pInt23.JPa - R1*pInt12.JVa;  //Equ 19 in main papaer
    M_k.block<3,3>(0,q) = -0.5*(pInt12.dT + pInt23.dT)*Eigen::Matrix3d::Identity();  //Equ 20 in main paper

    Eigen::Vector3d pi_k;
    pi_k = B*pInt23.GetDeltaPosition(bg, ba) - A*pInt12.GetDeltaPosition(bg, ba) + R1*pInt12.GetDeltaVelocity(bg, ba)
           + (T2.linear()-T1.linear())*Tcb.translation()/pInt12.dT
           - (T3.linear()-T2.linear())*Tcb.translation()/pInt23.dT;   //Equ 22 in main paper

    Eigen::Matrix3d Covariance;
    Covariance  = A*pInt12.C.block<3, 3>(6, 6)*A.transpose();
    Covariance += B*pInt23.C.block<3, 3>(6, 6)*B.transpose();
    Covariance += T1.linear()*pInt12.C.block<3, 3>(3, 3)*T1.linear().transpose();

    Eigen::Matrix3d Information = selfAdjointInverse(Covariance);
    //Eigen::Matrix3d Information = Eigen::Matrix3d::Identity();
    M +=  M_k.transpose()*Information*M_k;  //Equ 25 in main paper
    m += -2.*M_k.transpose()*Information*pi_k;  //Equ 25
    Q +=  pi_k.transpose()*Information*pi_k;  //Equ 25
  }

  //Solve
  Eigen::Matrix4d A = 2.*M.block<4, 4>(0, 0);

  //TODO check if A is invertible
  Eigen::Matrix3d A_ = A.block<3,3>(1,1);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> svdA_(A_, Eigen::EigenvaluesOnly);
  result.svA_ = svdA_.eigenvalues();
  result.detA_ = A_.determinant();

  Eigen::MatrixXd Bt = 2.*M.block<3, 4>(q, 0);
  Eigen::MatrixXd BtAi = Bt*A.inverse();

  Eigen::Matrix3d D = 2.*M.block<3, 3>(q, q);
  Eigen::Matrix3d S = D - BtAi*Bt.transpose();

  // TODO Check if S is invertible!
  //Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> svdS(S, Eigen::EigenvaluesOnly);
  //result.svS = svdS.eigenvalues();
  //result.detS = S.determinant();
  //LOG(INFO) << StringPrintf("det(S): %.16f", S.determinant());
  //LOG(INFO) << StringPrintf("eigenvalues(S): %.16f %.16f %.16f",c[0], svd.eigenvalues()[1], svd.eigenvalues()[2]);
  
  Eigen::Matrix3d Sa = S.determinant()*S.inverse();
  Eigen::Matrix3d U = S.trace()*Eigen::Matrix3d::Identity() - S;

  Eigen::Vector3d v1 = BtAi*m.head<q>();
  Eigen::Vector3d m2 = m.tail<3>();

  Eigen::Matrix3d X; Eigen::Vector3d Xm2;


}

