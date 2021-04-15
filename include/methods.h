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

  //This is the proposed analytical methods described in the paper
  //It uses Lagrange's multiplier
  
  //For detailed understanding read the appendix in the main paper and also
  //A. Censi, “An ICP variant using a point-to-line metric,” in IEEE Int.Conf. Robot. Autom. (ICRA), 2008.

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

  //Below we are solving the quadratic system with a quadratic constraint using Lagrange's multiplier
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

  X = I
  const double c4 = 16.*(v1.dot(  v1) - 2.*v1.dot( m2) + m2.dot( m2));

  X = U; Xm2 = X*m2;
  const double c3 = 16.*(v1.dot(X*v1) - 2.*v1.dot(Xm2) + m2.dot(Xm2));

  X = 2.*Sa + U*U; Xm2 = X*m2;
  const double c2 =  4.*(v1.dot(X*v1) - 2.*v1.dot(Xm2) + m2.dot(Xm2));

  X = Sa*U + U*Sa; Xm2 = X*m2;
  const double c1 =  2.*(v1.dot(X*v1) - 2.*v1.dot(Xm2) + m2.dot(Xm2));

  X = Sa*Sa; Xm2 = X*m2;
  const double c0 =     (v1.dot(X*v1) - 2.*v1.dot(Xm2) + m2.dot(Xm2));

  const double s00 = S(0, 0), s01 = S(0, 1), s02 = S(0, 2);
  const double s11 = S(1, 1), s12 = S(1, 2), s22 = S(2, 2);

  const double t1 = s00 + s11 + s22;
  const double t2 = s00*s11 + s00*s22 + s11*s22
                    - std::pow(s01, 2) - std::pow(s02, 2) - std::pow(s12, 2);

  const double t3 = s00*s11*s22 + 2.*s01*s02*s12
                    - s00*std::pow(s12, 2) - s11*std::pow(s02, 2) - s22*std::pow(s01, 2);

  Eigen::VectorXd coeffs(7);
  coeffs << 64.,
            64.*t1,
            16.*(std::pow(t1, 2) + 2.*t2),
            16.*(t1*t2 + t3),
            4.*(std::pow(t2, 2) + 2.*t1*t3),
            4.*t3*t2,
            std::pow(t3, 2);

  const double G2i = 1. / std::pow(IMU::GRAVITY_MAGNITUDE, 2);

  coeffs(2) -= c4*G2i;
  coeffs(3) -= c3*G2i;
  coeffs(4) -= c2*G2i;
  coeffs(5) -= c1*G2i;
  coeffs(6) -= c0*G2i;

  Eigen::VectorXd real,imag;
  if(!FindPolynomialRootsCompanionMatrix(coeffs,&real,&imag)){
    LOG(INFO) << "Failed to find the roots\n" 
              << StringPrintf("%.16f %.16f %.16f %.16f %.16f %.16f %.16f",
                               coeffs[0], coeffs[1], coeffs[2], coeffs[3],coeffs[4], coeffs[5], coeffs[6]);

    result.success = false;
    result.solve_ns = timer.ElapsedNanoSeconds();
    return;
  }

  Eigen::VectorXd lambdas = real_roots(real,imag);
  if(lambdas.size() == 0){
    LOG(INFO) << "No real roots found\n"
              << StringPrintf("%.16f %.16f %.16f %.16f %.16f %.16f %.16f",
                               coeffs[0], coeffs[1], coeffs[2], coeffs[3],coeffs[4], coeffs[5], coeffs[6]);

    result.success = false;
    result.solve_ns = timer.ElapsedNanoSeconds();
    return;
  }

  Eigen::MatrixXd W(n,n);
  W.setZero();
  W.block<3, 3>(q, q) = Eigen::Matrix3d::Identity();

  Eigen::VectorXd solution;
  double min_cost = std::numeric_limits<double>::max();
  for(Eigen::VectorXd::Index i=0; i<lambdas.size(); i++){
    const double lambda = lambdas(i);
    Eigen::FullPivLU<Eigen::MatrixXd> lu(2.*M + 2.*lambda*W);
    Eigen::VectorXd x_ = -lu.inverse()*m;

    double cost = x_.transpose()*M*x_;
    cost += m.transpose()*x_;
    cost += Q;

    if (cost < min_cost) {
      solution = x_;
      min_cost = cost;
    }
  }

  result.solve_ns = timer.ElapsedNanoSeconds();

  const double constraint = solution.transpose()*W*solution;
  if (solution[0] < 1e-3 || constraint < 0.|| std::abs(std::sqrt(constraint)-IMU::GRAVITY_MAGNITUDE)/IMU::GRAVITY_MAGNITUDE > 1e-3){
    LOG(INFO) << "Discarding bad solution...\n"
              << StringPrintf("scale: %.16f\n", solution[0])
              << StringPrintf("constraint: %.16f\n", constraint)
            << StringPrintf("constraint error: %.2f %",100.*std::abs(std::sqrt(constraint)-IMU::GRAVITY_MAGNITUDE)/IMU::GRAVITY_MAGNITUDE);
    
    result.success = false;
    return;
  }

  result.success = true;
  result.scale = solution[0];
  result.bias_a = solution.segment<3>(1);
  result.gravity = solution.segment<3>(4);
}

void iterative(const InputType& input, ResultType& result, const double initial_scale = 1.,
               const Eigen::Isometry3d& Tcb = Eigen::Isometry3d::Identity(),double* cost = nullptr,
               bool use_prior = true, double prior = 1e5){
  std::vector<double*> pointers;
  std::vector<double**> pointers2;

  //Global Parameters
  double* bg_ptr = new double[3];
  pointers.push_back(bg_ptr);
  Eigen::Map<Eigen::Vector3d> bg(bg_ptr);
  bg.setZero();

  double* ba_ptr = new double[3];
  pointers.push_back(ba_ptr);
  Eigen::Map<Eigen::Vector3d> ba(ba_ptr);
  ba.setZero();

  double* Rwg_ptr = new double[9];
  pointers.push_back(Rwg_ptr);
  Eigen::Map<Eigen::Matrix3d> Rwg(Rwg_ptr);
  Rwg.setIdentity();

  double* s_ptr = new double[1];
  pointers.push_back(s_ptr);
  s_ptr[0] = initial_scale;

  //Local Parameters (for each keyframe)
  double* v0_ptr = new double[3];
  pointers.push_back(v0_ptr);
  Eigen::Map<Eigen::Vector3d> v0(v0_ptr);
  v0.setZero();

  double** parameters = new double*[6];
  pointers2.push_back(parameters);
  parameters[0] = v0_ptr;
  parameters[1] = nullptr;
  parameters[2] = bg_ptr;
  parameters[3] = ba_ptr;
  parameters[4] = Rwg_ptr;
  parameters[5] = s_ptr;

  ceres::Problem problem;

  Eigen::Vector3d dirG;
  dirG.setZero();

  for(unsigned i=0; i<input.size(); i++){
    const Eigen::Isometry3d T1 = input[i].T1;
    const Eigen::Isometry3d T2 = input[i].T2;

    const std::shared_ptr<IMU::Preintegration> pInt = input[i].pInt;

    double* v1_ptr = parameters[0];
    Eigen::Map<Eigen::Vector3d> v1(v1_ptr);

    double* v2_ptr = new double[3];
    pointers.push_back(v2_ptr);
    Eigen::Map<Eigen::Vector3d> v2(v2_ptr);

    v2 = (T2.translation() - T1.translation())/pInt->dT;
    v1 = v2;

    parameters[1] = v2_ptr;

    //Rwg initialization
    dirG = -T1.linear()*pInt->GetUpdatedDeltaVelocity();

    ceres::CostFunction* cost_fucntion = new InertialCostFunction(pInt, T1.linear(),T1.translation(),
                                                                  T2.linear(), T2.translation(),Tcb);

    problem.AddResidualBlock(cost_fucntion,nullptr,parameters,6);

    double** parameters_ = new double*[6];
    pointers2.push_back(parameters_);
    parameters_[0] = parameters[1];
    parameters_[1] = nullptr;
    parameters_[2] = bg_ptr;
    parameters_[3] = ba_ptr;
    parameters_[4] = Rwg_ptr;
    parameters_[5] = s_ptr;

    parameters = parameters_;
  }

  if(use_prior){
    ceres::CostFunction* prior_cost_fucntion = new BiasPriorCostFunction(prior);
    problem.AddResidualBlock(prior_cost_function,nullptr,bg_ptr);
  }

  //Initialize Rwg estimate
  dirG = dirG.normalized();
  const Eigen::Vector3d gI = IMU::GRAVITY_VECTOR.normalized();
  const Eigen::Vector3d v = gI.cross(dirG);
  const double cos_theta = gI.dot(dirG);
  const double theta = std::acos(cos_theta);
  Rwg = ExpSO3(v*theta/v.norm());

  //Add Local parameterization
  GravityParameterization* gravity_local_parameterization = new GravityParameterization;
  problem.setParameterization(Rwg_ptr,gravity_local_parameterization);

  ScaleParameterization* scale_local_parameterization = new ScaleParameterization;
  problem.setParameterization(s_ptr,scale_local_parameterization);

  ceres::Solver::Options options;
  options.max_num_iterations = 200;
  ceres::Solver::Summary summary;

  Timer timer;
  timer.Start();

  ceres::Solve(options, &problem, &summary);

  timer.Pause();

  bool converged = (summary.termination_type == ceres::CONVERGENCE);
  if(converged){
    result.success = true;
    result.solve_ns = timer.ElapsedNanoSeconds();
    result.scale = s_ptr[0];
    result.bias_g = bg;
    result.bias_a = ba;
    result.gravity = Rwg*IMU::GRAVITY_VECTOR;

    if(cost){
      *cost = summary.final_cost;
    }
  }else{
    result.success = false;
    result.solve_ns = timer.ElapsedNanoSeconds();
  }

  //Free Memory
  for (double* ptr : pointers)
    delete[] ptr;
  for (double** ptr : pointers2)
    delete[] ptr;

}

