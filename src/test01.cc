// 10KFs at 4Hz

#define PROGRAM_NAME \
  "test01"

#define FLAG_CASES \
  FLAG_CASE(string, logs_dir, "./logs/", "Logs save directory")

#define ARGS_CASES \
  ARG_CASE(dataset_dir)

// STL
#include <algorithm>
#include <cmath>


#include <iterator>
#include <string>
#include <vector>

//Boost
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>

//Eigen
#include <Eigen/Core>
#include <Eigen/Geometry>

//Glog
#include <glog/logging.h>

#include "imu_preintegration.h"
#include "io.h"
#include "methods.h"

#include "util/args.h"
#include "util/csv.h"

namespace fs = boost::filesystem;

using GroundTruth = std::vector<io::state_t>;
using ImuData = io::ImuData;

//IMU parameters
//Euroc dataset

const double rate = 200;
const double dt = 1./rate;
const double ng = 1.7e-4;
const double na = 2.e-3;

struct evaluation_t{
  evaluation_t(const std::uint64_t solve_time, const std::uint64_t initialization_time,const std::uint64_t timestamp,
               const double scale_error, const double gyro_bias_error, const double acc_bias_error, const double gravity_error) : 
    solve_time(solve_time), initialization_time(initialization_time), timestamp(timestamp), scale_error(scale_error), gyro_bias_error(gyro_bias_error),
    acc_bias_error(acc_bias_error), gravity_error(gravity_error) {}

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  std::uint64_t solve_time;  //nanoseconds
  std::uint64_t initialization_time;  //nanoseconds
  std::uint64_t timestamp;  //nanoseconds
  const double scale_error;  //percent
  const double gyro_bias_error;  //percent
  const double acc_bias_error;  //percent
  const double gravity_error;  //degrees

  //Analysis
  Eigen::Vector4d svA;
  double detA;
  Eigen::Vector3d svS;
  double detS;
}; //evaluation_T

void ValidateArgs(){
  CHECK(fs::is_directory(ARGS_dataset_dir));
}

void ValidateFlags(){
  fs::create_directories(FLAGS_logs_dir);
}

GroundTruth::const_iterator next(GroundTruth::const_iterator i, GroundTruth::const_iterator j, io::timestamp_t dt){
  if(i == j){
    LOG(WARNING) << "Already at the end...";
    return i;
  }

  io::timestamp_t t = i->timestamp + dt;

  GroundTruth::const_iterator it = std::lower_bound(i,j,t,[](const io::state_t& lhs, const io::timestamp_t rhs){return lhs.timestamp < rhs;});

  if(it == i)
    return i;
  if(it == j)
    return j;

  GroundTruth::const_iterator it_ = std::next(it,-1);

  if((it->timestamp - t) > (t - it_->timestamp))
    return it_;
  else
    return it;
}

ImuData::const_iterator start_imu(ImuData::const_iterator i, ImuData::const_iterator j, io::timestamp_j t){

  ImuData::const_iterator = std::upper_bound(i,j,t,[](const io::timestamp_t lhs, const io::imu_data_t& rhs){return lhs < rhs.timestamp;});

  if (it == i) return i;
  if (it == j) return j;
  ImuData::const_iterator it_ = std::next(it, -1);

  if ((it->timestamp - t) > (t - it_->timestamp))
    return it_;
  else
    return it;
}

GroundTruth::const_iterator start(const GroundTruth& trajectory, const io::ImuData& imu_data, const unsigned nframes){
  GroundTruth::const_iterator i = trajectory.cbegin();
  GroundTruth::const_iterator i_ = i;

  while(i != trajectory.cend()){
    Eigen::Vector3d angA;
    avgA.setZero();

    io::ImuData::const_iterator it = imu_data.cbegin();
    for (unsigned n = 0; n < nframes; ++n){
      it = start_imu(it, imu_data.cend(), i->timestamp);
      CHECK(it != imu_data.cend());

      Groundtruth::const_iterator j = next(i, trajectory.cend(), 250000000); // 4 Hz
      CHECK(j != trajectory.cend());

      std::shared_ptr<IMU::Preintegrated> pInt = std::make_shared<IMU::Preintegrated>(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
      while(it!= imu_data.cend() && std::llabs(it->timestamp - j->timestamp) > 2500000){
        const Eigen::Vector3d w(it->w_x, it->w_y, it->w_z);
        const Eigen::Vector3d a(it->a_x, it->a_y, it->a_z);
        pInt->IntegrateNewMeasurement(w, a, dt);
        std::advance(it, 1);
      }
      CHECK(it != imu_data.cend());
      avgA += pInt->dV/pInt->dT;
      i=j;
    }

    avgA /= static_cast<double>(nframes);
    const double avgA_error = std::abs(avgA.norm() - IMU::GRAVITY_MAGNITUDE) / IMU::GRAVITY_MAGNITUDE;
    LOG(INFO) << "Average acceleration: " << 100.*avgA_error;
    if (avgA_error > 5e-3) 
      break;

    i = next(i, trajectory.cend(), 500000000); // 0.5s
    i_ = i;
  }
  CHECK(i != trajectory.cend());

  return i_;
}

void save(std::vector<evaluation_t>& data, std::string& save_path){
  Eigen::MatrixXd m(data.size(),7);

  for(unsigned i=0; i<data.size(); i++){
    Eigen::RowVectorXd row(7);
    row << data[i].solve_time,data[i].initialization_time, data[i].timestamp,
           data[i].scale_error, data[i].gyro_bias_error, data[i].acc_bias_error, data[i].gravity_error;

    m.row(i) = row;
  }

  csv::write(m,save_path);
}

void saveExt(std::vector<evaluation_t>& data, std::string& save_path){
  Eigen::MatrixXd m(data.size(), 8);

  for(unsigned i=0; i<data.size(); i++){
    Eigen::RowVectorXd row(8);
    row << data[i].solve_time, data[i].initialization_time, data[i].timestamp,
        //4
           data[i].scale_error, data[i].gyro_bias_error, data[i].acc_bias_error, data[i].gravity_error,
        //8
           //data[i].svA(0), data[i].svA(1), data[i].svA(2), data[i].svA(3), data[i].detA,
        //13
           //data[i].svS(0), data[i].svS(1), data[i].svS(2), data[i].detS;
    m.row(i) = row;
  }

  csv::write(m,save_path);
}

void run(const fs::path& sequence_path){
  std::string sequence_name = sequence_path.filename().string();
  if(sequence_name == "."){
    sequence_name = sequence_path.parent_path().filename().string();
  }

  LOG(INFO) << "Running Experiment" << sequence_name;

  fs::path trajectory_path = sequence_path / "state_groundtruth_estimate0" / "data.csv";
  CHECK(fs::is_regular_file(trajectory_path)) << "Path not found: " << trajectory_path.string();

  fs::path data_path = sequence_path/"imu0"/"data.csv";
  CHECK(fs::is_regular_file(data_path)) << "Path not found" << data_path.string();

  GroundTruth trajectory = io::read_file<GroundTruth::value_type>(trajectory_path.string());
  io::ImuData imu_data = io::read_file<io::ImuData::value_type>(data_path.string());

  std::vector<unsigned> possible_nframes = {5,10,20,50,75};

  for(unsigned nframes : possible_nframes){
    GroundTruth::const_iterator i = start(trajectory,imu_data,nframes);
    LOG(info) << "Starting at: " << i->timestamp;
    LOG(info) << StringPrintf("With %d frames",nframes);

    std::vector<evaluation_t> proposed_evaluation;
    std::vector<evaluation_t> iterative_evaluation;

    unsigned count = 0;
    std::uint64_t imu_integration = 0;
    std::uint64_t skipped = 0;

    GroundTruth::const_iterator i_ = i;
    while(i != trajectory.cend()){
      InputType input;

      Eigen::Vector3d avgBg = Eigen::Vector3d(i->bw_x,i->bw_y,i->bw_z);
      Eigen::Vector3d avgBa = Eigen::Vector3d(i->ba_x,i->ba_y,i->ba_z);

      Eigen::Vector3d avgA;
      avgA.setZero();

      io::ImuData::const_iterator it = imu_data.cbegin();
      for(unsigned n=0; n<nframes; n++){

        it = start_imu(it, imu_data.cend(), i->timestamp);
        if (it == imu_data.cend()) {
          LOG(WARNING) << "Couldn't find IMU measurement at " << i->timestamp;
          break;
        }

        Groundtruth::const_iterator j = next(i, trajectory.cend(), 250000000); // 4 Hz
        if (j == trajectory.cend()) {
          LOG(WARNING) << "Couldn't find next frame for " << i->timestamp;
          break;
        }

        Timer timer;
        timer.start();

        std::shared_ptr<IMU::Preintegration> pInt = std::make_shared<IMU::Preintegration>(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
        while(it != imu_data.cend() && std::llabs(it->timestamp - j->timestamp) > 2500000){
          const Eigen::Vector3d w(it->w_x, it->w_y, it->w_z);
          const Eigen::Vector3d a(it->a_x, it->a_y, it->a_z);
          pInt->IntegrateNewMeasurement(w, a, dt);
          std::advance(it, 1);
        }  //while loop line 254

        if(it == imu_data.cend()){
          LOG(WARNING) << "IMU stream ended!";
          break;
        }

        imu_integration += timer.ElapsedNanoSeconds();
        count++;

        avgBg += Eigen::Vector3d(j->bw_x, j->bw_y, j->bw_z);
        avgBa += Eigen::Vector3d(j->ba_x, j->ba_y, j->ba_z);

        avgA += pInt->dV/pInt->dT;
        input.emplace_back(i->pose, i->timestamp, j->pose, j->timestamp, pInt);

        i=j;

      } //for loop line 236

      if (input.size() < nframes) {
        LOG(INFO) << StringPrintf("I don't have %d frames. I think dataset ended...", nframes);
        break;
      }

      avgBg /= static_cast<double>(nframes+1);
      avgBa /= static_cast<double>(nframes+1);

      avgA /= static_cast<double>(nframes);
      const double avgA_error = std::abs(avgA.norm() - IMU::GRAVITY_MAGNITUDE) / IMU::GRAVITY_MAGNITUDE;

      if(avgA_error > 5e-3){
        std::uint64_t timestamp = input[0].t1;
        std::uint64_t initialization_time = i->timestamp - i_->timestamp;

        {
          ResultType gyroscope_result;
          proposed_gyroscope(input, gyroscope_result);

          ResultType accelerometer_result;
          proposed_accelerometer(input, accelerometer_result, gyroscope_result.bias_g);

          ResultType proposed_result;
          proposed_result.success = gyroscope_result.success & accelerometer_result.success;
          proposed_result.solve_ns = gyroscope_result.solve_ns + accelerometer_result.solve_ns;
          proposed_result.scale = accelerometer_result.scale;
          proposed_result.bias_g = gyroscope_result.bias_g;
          proposed_result.bias_a = accelerometer_result.bias_a;
          proposed_result.gravity = accelerometer_result.gravity;

          if(proposed_result.success){
            const double scale_error = 100.*std::abs(proposed_result.scale - 1.);
            const double gyro_bias_error = 100.*std::abs(proposed_result.bias_g.norm() - avgBg.norm()) / avgBg.norm();
            const double acc_bias_error = 100.*std::abs(proposed_result.bias_a.norm() - avgBa.norm()) / avgBa.norm();
            const double gravity_error = 180.*std::acos(proposed_result.gravity.normalized().dot(IMU::GRAVITY_VECTOR.normalized()))/EIGEN_PI;
            evaluation_t eval_data(proposed_result.solve_ns, initialization_time, timestamp,
                                   scale_error, gyro_bias_error, acc_bias_error, gravity_error);
            proposed_evaluation.push_back(eval_data);
          }else{
            LOG(ERROR) << "Proposed method failed at " << timestamp;
          }
        }

        {
          ResultType iterative_result;
          iterative(input,iterative_Result);

          if(iterative_result.success){
            const double scale_error = 100.*std::abs(iterative_result.scale - 1.);
            const double gyro_bias_error = 100.*std::abs(iterative_result.bias_g.norm() - avgBg.norm()) / avgBg.norm();
            const double acc_bias_error = 100.*std::abs(iterative_result.bias_a.norm() - avgBa.norm()) / avgBa.norm();
            const double gravity_error = 180.*std::acos(iterative_result.gravity.normalized().dot(IMU::GRAVITY_VECTOR.normalized()))/EIGEN_PI;
            iterative_evaluation.emplace_back(iterative_result.solve_ns, initialization_time, timestamp,
                                              scale_error, gyro_bias_error, acc_bias_error, gravity_error);
          }else{
            LOG(ERROR) << "Iterative method failed at" << timestamp; 
          }
        }

        i = next(i_, trajectory.cend(), 500000000);
        i_ = i;
        skipped = 0;
      }else{
        skipped += 500000000;
        i = next(i_, trajectory.cend(), skipped); // 0.5s
      }
    } //while loop line 226

    std::string proposed_file = StringPrintf("%s_%d_ours.csv", sequence_name.c_str(), nframes);
    LOG(INFO) << "Saving evaluation data into " << proposed_file;
    save(proposed_evaluation, proposed_file);

    std::string iterative_file = StringPrintf("%s_%d_iterative.csv", sequence_name.c_str(), nframes);
    LOG(INFO) << "Saving evaluation data into " << iterative_file;
    save(iterative_evaluation, iterative_file);
  }  //for loop line 213

  LOG(INFO) << "done." << std::endl;

}

int main(int argc, char* argv[]){
  //Handle help flag
  if(args::HelpRequired(argc,argv)){
    argc::ShowHelp();
    return 0;
  }

  //Parse input flag
  args::ParseCommandLineNonHelpFlags(&argc,&argv,true);

  FLAGS_log_dir = FLAGS_logs_dir;
  FLAGS_stderrthreshold = 0;
  google::InitGoogleLogging(argv[0]);

  //check number of args
  if(argc-1 != args::NumArgs()){
    args::ShowHelp();
    return -1;
  }

  //Parse input args
  args::ParseCommandLineArgs(argc, argv);

  //Validate input arguments
  ValidateFlags();
  ValidateArgs();

  IMU::Sigma.block<3, 3>(0, 0) = rate*ng*ng * Eigen::Matrix3d::Identity();
  IMU::Sigma.block<3, 3>(3, 3) = rate*na*na * Eigen::Matrix3d::Identity();

  run(ARGS_dataset_dir);

  return 0;

}
