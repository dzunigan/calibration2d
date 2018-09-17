// Calibration

// STL
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

// Boost
#include <boost/filesystem.hpp>
#include <boost/math/special_functions/sign.hpp>

// Glog
#include <glog/logging.h>

#include "loransac.h"
#include "estimators.hpp"

#include "refinement.hpp"

#include "util/alignment.h"
#include "util/csv.hpp"

// TODO: assuming the last sensor is the only nonmetrically accurate (for additional constraints)

constexpr size_t MIN_INLIERS = 5;
#define NUM_ITERATIONS 1

using Calibration_RANSAC = colmap::LORANSAC<CalibrationEstimator2, CalibrationEstimator2>;
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_CUSTOM(Calibration_RANSAC::Report)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_CUSTOM(Sensor)

#ifndef PROGRAM_NAME
#define PROGRAM_NAME \
    "calibrate"
#endif

inline void PrintHelp(std::ostream& out) {
    out << "Usage: " << PROGRAM_NAME;
    out << " <motions1> <motions2> ...";
    out << std::endl;
}

inline void PrintCSV(std::ostream& out, const Eigen::Ref<const Eigen::Vector4d>& vec, int precision = 6) {
    out << std::fixed << std::setprecision(precision)
        << vec(0) << ", " << vec(1) << ", " << vec(2) << ", " << vec(3);
}

int main(int argc, char* argv[]) {

    // Initialize Google's logging library
    google::InitGoogleLogging(argv[0]);
    FLAGS_logtostderr = true;

    // Check at runtime that size_t can contain all sensor ids without information loss
    CHECK_EQ(static_cast<size_t>(kMaxNumSensors), kMaxNumSensors);

    // At least two sensors required
    if (argc < 3) {
        PrintHelp(std::cerr);
        return -1;
    }

    // Check input args
    for (int i = 1; i < argc; ++i)
        CHECK(boost::filesystem::is_regular_file(argv[i])) << "Invalid input file:" << std::endl
                                                           << argv[i];

    // Load sensor observations
    std::vector<Sensor> sensors(argc-1);
    CHECK_GE(sensors.size(), 2) << "At least 2 sensors required to perform extrinsic calibration!";
    for (int i = 1; i < argc; ++i) {
        Eigen::MatrixXd data = csv::read<double>(argv[i], ' ');

        Sensor& sensor = sensors[i-1];
        sensor.SetSensorId(i-1);

        for (int k = 0; k < data.rows(); ++k)
            sensor.Observations().push_back(data.row(k));
    }

    // Set reference sensor
    const Sensor& reference = sensors.front();
    size_t N = reference.Observations().size();

    // All sensor must have the same number of observations
    for (const Sensor& sensor : sensors)
        CHECK_EQ(sensor.Observations().size(), N);

//    const size_t NUM_CALIBRATIONS = sensors.size() - 1;
//    Eigen::MatrixXd data(NUM_ITERATIONS, 3*4*NUM_CALIBRATIONS);
//    data.fill(std::numeric_limits<double>::quiet_NaN());

    // Outputs
    std::vector<Eigen::Vector4d> closed_form;
    std::vector<Eigen::Vector4d> iterative;
    std::vector<Eigen::Vector4d> robust;

    // Auto-initialization
    // -------------------
    ObservationDatabase observations; // RANSAC inliers

    colmap::RANSACOptions ransac_options;
    ransac_options.max_error = 0.1;
    ransac_options.min_inlier_ratio = 0.25;
    ransac_options.confidence = 0.99999;
    ransac_options.min_num_trials = 100;
    ransac_options.max_num_trials = 1000;

    std::vector<Calibration_RANSAC::Report> reports(sensors.size() - 1);
    for (size_t i = 1; i < sensors.size(); ++i) {
        Sensor& sensor = sensors[i];
        Calibration_RANSAC::Report& ransac_report = reports[i-1];

        Calibration_RANSAC calibration_ransac(ransac_options);
        ransac_report = calibration_ransac.Estimate(reference.Observations(), sensor.Observations());

        CHECK_EQ(sensor.Observations().size(), ransac_report.inlier_mask.size());
        for (size_t k = 0; k < sensor.Observations().size(); ++k) {
            if (ransac_report.inlier_mask.at(k))
                observations.AddIndexPair(reference.SensorId(), sensor.SensorId(), k, k);
        }

        // Check inlier set
        CHECK(ransac_report.inlier_mask.size() > ransac_options.min_inlier_ratio * sensor.Observations().size()) << "Not enough inliers!";
        sensor.Calibration() = ransac_report.model;
        closed_form.push_back(ransac_report.model);

        if (i != 1) std::cout << ", ";
        // std::cout << ransac_report.model.transpose();
        PrintCSV(std::cout, ransac_report.model);
    }

    // Additional constraints
    std::unordered_set<index_pair_t> constraints;
    for (size_t i = 1; i < sensors.size(); ++i) {
        const sensor_t id_i = sensors[i].SensorId();
        const Calibration_RANSAC::Report& ransac_report_i = reports[i-1];
        for (size_t j = i+1; j < sensors.size(); ++j) {
            const size_t id_j = sensors[j].SensorId();
            const Calibration_RANSAC::Report& ransac_report_j = reports[j-1];

            CHECK_EQ(ransac_report_i.inlier_mask.size(), ransac_report_j.inlier_mask.size());
            for (size_t k = 0; k < N; ++k) {
                if (ransac_report_i.inlier_mask.at(k) && ransac_report_j.inlier_mask.at(k))
                    observations.AddIndexPair(id_i, id_j, k, k); // Intersection of inliers
            }

            if (observations.Pairs(id_i, id_j).size() > MIN_INLIERS)
                 constraints.insert(index_pair_t(id_i, id_j)); // TODO: In metric space (in first place the sensor providing metrically accurate incremental poses)
        }
    }

    // Refinement options
    BatchCalibrationOptions options;
    options.solver_options.minimizer_progress_to_stdout = false;
    options.print_summary = false;

    // Refinement configuration
    BatchCalibrationConfig config;
    for (const Sensor& sensor : sensors)
        config.AddSensor(sensor);

    CHECK(config.HasSensor(reference.SensorId()));
    config.SetReferenceSensor(reference.SensorId());
    config.ObservationDatabase() = observations;
    config.AdditionalConstraints() = constraints;

    // Iterative
    // ---------
    options.use_additional_constraints = false;
    options.loss_function_type = BatchCalibrationOptions::LossFunctionType::TRIVIAL;

    // Run iterative optimization
    BatchCalibration calibration_iterative(options, config);

    if (calibration_iterative.Solve()) {
        for (size_t i = 1; i < sensors.size(); ++i) {
            const Sensor& sensor = sensors[i];
            Eigen::Vector4d params = calibration_iterative.Paramerameters(sensor.SensorId());
            params *= boost::math::sign(params(3));

            iterative.push_back(params);
            std::cout << ", "; // << params;
            PrintCSV(std::cout, params);
        }
    } else {
        Eigen::Vector4d vnan;
        vnan.fill(std::numeric_limits<double>::quiet_NaN());
        for (size_t i = 1; i < sensors.size(); ++i) {
            std::cout << ", "; // << vnan.transpose() << std::endl;
            PrintCSV(std::cout, vnan);
        }
    }

    // Robust
    // ------
    options.use_additional_constraints = true;
    options.loss_function_type = BatchCalibrationOptions::LossFunctionType::CAUCHY;
    options.loss_function_scale = 0.05;

    // Run robust optimization
    BatchCalibration calibration_robust(options, config);

    if (calibration_robust.Solve()) {
        for (size_t i = 1; i < sensors.size(); ++i) {
            const Sensor& sensor = sensors[i];
            Eigen::Vector4d params = calibration_robust.Paramerameters(sensor.SensorId());
            params *= boost::math::sign(params(3));

            robust.push_back(params);
            std::cout << ", "; // << params;
            PrintCSV(std::cout, params);
        }
    } else {
        Eigen::Vector4d vnan;
        vnan.fill(std::numeric_limits<double>::quiet_NaN());
        for (size_t i = 1; i < sensors.size(); ++i) {
            std::cout << ", "; // << vnan.transpose() << std::endl;
            PrintCSV(std::cout, vnan);
        }
    }

    // Save closed_form, iterative, robust estimations
//    for (size_t k = 0; k < NUM_CALIBRATIONS; ++k) {
//        if (k < closed_form.size())
//            data.block<1, 4>(asdf, 4*k) = closed_form[k];
//        if (k < iterative.size())
//            data.block<1, 4>(asdf, 4*NUM_CALIBRATIONS + 4*k) = iterative[k];
//        if (k < robust.size())
//            data.block<1, 4>(asdf, 8*NUM_CALIBRATIONS + 4*k) = robust[k];
//    }

    std::cout << std::endl;

    // Output calibrations
//    CHECK(csv::write(data, "data.csv"));

    return 0;
}
