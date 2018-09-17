#include "estimators.hpp"

// Boost
#include <boost/math/special_functions/sign.hpp>

// Eigen
#include <Eigen/Dense>

#include "util/logging.h"
#include "util/math.h"

#include "sim2.hpp"

/*
std::vector<RotationAndScaleEstimator::M_t>
RotationAndScaleEstimator::Estimate(
        const std::vector<X_t>& s1, const std::vector<Y_t>& s2) {
    CHECK_GE(s1.size(), 1);
    CHECK_GE(s2.size(), 1); // Not really needed
    CHECK_EQ(s1.size(), s2.size());

    Eigen::Matrix3d M = Eigen::Matrix3d::Zero();

    const std::size_t n = s1.size();
    for (std::size_t k = 0; k < n; ++k) {
        const X_t &x = s1[k];
        const Y_t &y = s2[k];

        Eigen::Matrix<double, 2, 3> Q;
        Q << -x(0), y(0), -y(1),
             -x(1), y(1),  y(0);

        M += Q.transpose()*Q;
    }

    std::cout << "M:" << std::endl;
    std::cout << M << std::endl;

    const double m00 = M(0, 0),
            m01 = 0.5*(M(0, 1) + M(1, 0)),
            m02 = 0.5*(M(0, 2) + M(2, 0)),
            m11 = M(1, 1);

    const double m01_2 = m01*m01,
            m02_2 = m02*m02;

    const double a = m00,
            b = - m01_2 - m02_2 + 2.0*m00*m11,
            c = -m11*(m01_2 + m02_2 - m00*m11);

    std::cout << "a: " << a << std::endl;
    std::cout << "b: " << b << std::endl;
    std::cout << "c: " << c << std::endl;

    double x0, x1;
    const int real_solutions = solveQuadratic(a, b, c, x0, x1);
    std::cout << "Real solutions: " << real_solutions << std::endl;

    std::vector<M_t> models;

    if (real_solutions == 0)
        return models;

    Eigen::Matrix3d W = Eigen::Matrix3d::Zero();
    W(1, 1) = 1; W(2, 2) = 1;

    std::vector<double> lambdas;
    lambdas.push_back(x0);
    if (real_solutions > 1)
        lambdas.push_back(x1);

    M_t solution; // Unique solution
    bool solution_found = false;
    double solution_cost = std::numeric_limits<double>::max();

    std::vector<double> residuals;
    for (const double lambda : lambdas) {
        std::cout << "lambda: " << lambda << std::endl;
        Eigen::MatrixXd x;
        std::cout << "M + lambda*W" << std::endl;
        std::cout << M + lambda*W << std::endl;
        int kernel_size = solveNullspace(M + lambda*W, x);
        std::cout << "Kernel size: " << kernel_size << std::endl;
        if (kernel_size != 1) continue;

        x *= boost::math::sign(x(0, 0)) / x.block<2, 1>(1, 0).norm(); // Normalize solution

        M_t candidate_solution;
        candidate_solution << x(0), std::atan2(x(2), x(1));
        Residuals(s1, s2, candidate_solution, &residuals);

        std::cout << "candidate_solution: " << candidate_solution.transpose() << std::endl;

        double cost = std::accumulate(residuals.begin(), residuals.end(), 0.0);
        std::cout << "cost: " << cost << std::endl;
        if (cost < solution_cost) {
            solution = candidate_solution;
            solution_found = true;
            solution_cost = cost;
        }
    }

    if (solution_found)
        models.push_back(solution);

    return models;
}

void RotationAndScaleEstimator::Residuals(
        const std::vector<X_t>& s1, const std::vector<Y_t>& s2, const M_t& c,
        std::vector<double>* residuals) {
    CHECK_EQ(s1.size(), s2.size());

    const std::size_t n = s1.size();
    residuals->resize(n);

    const double s = c(0);
    const double ct = std::cos(c(1));
    const double st = std::sin(c(1));

    for (std::size_t k = 0; k < n; ++k) {
        const X_t &x = s1[k];
        const Y_t &y = s2[k];

        const double ex = ct*y(0) - st*y(1) - s*x(0);
        const double ey = st*y(0) + ct*y(1) - s*x(1);

        (*residuals)[k] = ex*ex + ey*ey;
    }
}
*/

std::vector<CalibrationEstimator::M_t>
CalibrationEstimator::Estimate(
        const std::vector<X_t>& x, const std::vector<Y_t>& y) {
    CHECK_GE(x.size(), kMinNumSamples);
    CHECK_GE(y.size(), kMinNumSamples); // Not really needed
    CHECK_EQ(x.size(), y.size());

    Eigen::Matrix5d M = Eigen::Matrix5d::Zero();

    const std::size_t n = x.size();
    for (std::size_t k = 0; k < n; ++k) {
        const X_t &x_k = x[k];
        const Y_t &y_k = y[k];

        double cx = std::cos(x_k(2)),
               sx = std::sin(x_k(2));

        Eigen::Matrix<double, 2, 5> Q;
        Q << -x_k(0), (1-cx),     sx, y_k(0), -y_k(1),
             -x_k(1),    -sx, (1-cx), y_k(1),  y_k(0);

        M += Q.transpose()*Q;
    }

//    std::cout << "M:" << std::endl;
//    std::cout << M << std::endl;

    const double m00 = M(0, 0),
            m01 = 0.5*(M(0, 1) + M(1, 0)),
            m02 = 0.5*(M(0, 2) + M(2, 0)),
            m03 = 0.5*(M(0, 3) + M(3, 0)),
            m04 = 0.5*(M(0, 4) + M(4, 0)),
            m11 = M(1, 1),
            m13 = 0.25*(M(1, 3) + M(2, 4) + M(3, 1) + M(4, 2)),
            m23 = 0.25*(M(2, 3) - M(1, 4) + M(3, 2) - M(4, 1)),
            m33 = M(3, 3);

    using std::pow;
    const double a = -m11*(pow(m01,2) + pow(m02,2) - m00*m11),
            b = -2.0*m33*pow(m01,2)*m11 + pow(m01,2)*pow(m13,2) + pow(m01,2)*pow(m23,2) + 2.0*m01*m03*m11*m13 - 2.0*m01*m04*m11*m23 - 2.0*m33*pow(m02,2)*m11 + pow(m02,2)*pow(m13,2) + pow(m02,2)*pow(m23,2) + 2.0*m02*m03*m11*m23 + 2.0*m02*m04*m11*m13 - pow(m03,2)*pow(m11,2) - pow(m04,2)*pow(m11,2) + 2.0*m00*m33*pow(m11,2) - 2.0*m00*m11*pow(m13,2) - 2.0*m00*m11*pow(m23,2),
            c = (pow(m13,2) + pow(m23,2) - m11*m33)*(m33*pow(m01,2) - 2.0*m01*m03*m13 + 2.0*m01*m04*m23 + m33*pow(m02,2) - 2.0*m02*m03*m23 - 2*m02*m04*m13 + m11*pow(m03,2) + m11*pow(m04,2) + m00*pow(m13,2) + m00*pow(m23,2) - m00*m11*m33);

//    std::cout << "a: " << a << std::endl;
//    std::cout << "b: " << b << std::endl;
//    std::cout << "c: " << c << std::endl;

    double x0, x1;
    const int real_solutions = solveQuadratic(a, b, c, x0, x1);

    std::vector<M_t> models;

    if (real_solutions == 0)
        return models;

    Eigen::Matrix5d W = Eigen::Matrix5d::Zero();
    W(3, 3) = 1; W(4, 4) = 1;
//    std::cout << "W:" << std::endl;
//    std::cout << W << std::endl;

    std::vector<double> lambdas;
    lambdas.push_back(x0);
    if (real_solutions > 1)
        lambdas.push_back(x1);

    M_t solution; // Unique solution
    bool solution_found = false;
    double solution_cost = std::numeric_limits<double>::max();

    std::vector<double> residuals;
    for (const double lambda : lambdas) {
//        std::cout << "lambda: " << lambda << std::endl;
        Eigen::MatrixXd phi;
//        std::cout << "M + lambda*W" << std::endl;
//        std::cout << M + lambda*W << std::endl;
        int kernel_size = solveNullspace(M + lambda*W, phi);
//        std::cout << "Kernel size: " << kernel_size << std::endl;
        if (kernel_size != 1) continue;

        phi *= boost::math::sign(phi(0, 0)) / phi.block<2, 1>(3, 0).norm(); // Normalize solution

        M_t candidate_solution;
        candidate_solution << phi(1), phi(2), std::atan2(phi(4), phi(3)), 1.0/phi(0);
        Residuals(x, y, candidate_solution, &residuals);

        double cost = std::accumulate(residuals.begin(), residuals.end(), 0.0);
//        std::cout << "cost: " << cost << std::endl;
        if (cost < solution_cost) {
            solution = candidate_solution;
            solution_found = true;
            solution_cost = cost;
        }
    }

    if (solution_found)
        models.push_back(solution);

    return models;
}

void CalibrationEstimator::Residuals(const std::vector<X_t>& x, const std::vector<Y_t>& y, const M_t& m,
        std::vector<double>* residuals) {
    CHECK_EQ(x.size(), y.size());
    CHECK_GT(m(3), 0);

    const std::size_t n = x.size();
    residuals->resize(n);

//    const double ct = std::cos(m(2));
//    const double st = std::sin(m(2));

//    const double s = m(3);

//    const double inv_x = -s*m(0)*ct - s*m(1)*st;
//    const double inv_y =  s*m(0)*st - s*m(1)*ct;

    for (std::size_t k = 0; k < n; ++k) {
        const X_t &x_k = x[k];
        const Y_t &y_k = y[k];

//        const double cy = std::cos(y_k(2));
//        const double sy = std::sin(y_k(2));

        Eigen::Vector4d y_k_(y_k(0), y_k(1), y_k(2), 1.0);
        Eigen::Vector4d e = PoseComposition<double>(m, PoseComposition<double>(y_k_, PoseInverse<double>(m)));
        const double ex = x_k(0) - e(0);
        const double ey = x_k(1) - e(1);
//        const double ex = s*y_k(0) + inv_x*cy - inv_y*sy - inv_x - x_k(0)*ct - x_k(1)*st;
//        const double ey = s*y_k(1) + inv_x*sy + inv_y*cy - inv_y + x_k(0)*st - x_k(1)*ct;
//        const double ex = m_x + y_k(0)*ct - y_k(1)*st - inv_s*x_k(0) - m_x*cx + m_y*sx;
//        const double ey = m_y + y_k(0)*st + y_k(1)*ct - inv_s*x_k(1) - m_x*sx - m_y*cx;

        (*residuals)[k] = ex*ex + ey*ey;
    }
}

std::vector<CalibrationEstimator2::M_t>
CalibrationEstimator2::Estimate(
        const std::vector<X_t>& x, const std::vector<Y_t>& y) {
    CHECK_GE(x.size(), kMinNumSamples);
    CHECK_GE(y.size(), kMinNumSamples); // Not really needed
    CHECK_EQ(x.size(), y.size());

    const std::size_t n = x.size();

    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(2 * n, 4);
    Eigen::MatrixXd w = Eigen::MatrixXd::Zero(2 * n, 1);
    for (std::size_t k = 0; k < n; ++k) {
        const X_t &x_k = x[k];
        const Y_t &y_k = y[k];

        Eigen::Matrix2d J;
        J = RotationMatrix2D(x_k(2)) - Eigen::Matrix2d::Identity();

        Eigen::Vector2d pi = y_k.head<2>();

        Eigen::Matrix2d K;
        K << pi(0), -pi(1), pi(1), pi(0);

        G.block<2,2>(k * 2, 0) = J;
        G.block<2,2>(k * 2, 2) = K;

        w.block<2,1>(k * 2, 0) = x_k.head<2>();
    }

    Eigen::MatrixXd m(4, 1);
    m = G.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(w);

    const double tx = -m(0);
    const double ty = -m(1);
    const double alpha = std::atan2(m(3), m(2));
    const double scale = m.block<2,1>(2, 0).norm();

    std::vector<M_t> models;

    M_t solution;
    solution << tx/scale, ty/scale, alpha, scale;

    models.push_back(solution);

    return models;
}

void CalibrationEstimator2::Residuals(const std::vector<X_t>& x, const std::vector<Y_t>& y, const M_t& m,
        std::vector<double>* residuals) {
    CHECK_EQ(x.size(), y.size());
    CHECK_GT(m(3), 0);

    const std::size_t n = x.size();
    residuals->resize(n);

//    const double ct = std::cos(m(2));
//    const double st = std::sin(m(2));

//    const double s = m(3);

//    const double inv_x = -s*m(0)*ct - s*m(1)*st;
//    const double inv_y =  s*m(0)*st - s*m(1)*ct;

    for (std::size_t k = 0; k < n; ++k) {
        const X_t &x_k = x[k];
        const Y_t &y_k = y[k];

//        const double cy = std::cos(y_k(2));
//        const double sy = std::sin(y_k(2));

        Eigen::Vector4d y_k_(y_k(0), y_k(1), y_k(2), 1.0);
        Eigen::Vector4d e = PoseComposition<double>(m, PoseComposition<double>(y_k_, PoseInverse<double>(m)));
        const double ex = x_k(0) - e(0);
        const double ey = x_k(1) - e(1);
//        const double ex = s*y_k(0) + inv_x*cy - inv_y*sy - inv_x - x_k(0)*ct - x_k(1)*st;
//        const double ey = s*y_k(1) + inv_x*sy + inv_y*cy - inv_y + x_k(0)*st - x_k(1)*ct;
//        const double ex = m_x + y_k(0)*ct - y_k(1)*st - inv_s*x_k(0) - m_x*cx + m_y*sx;
//        const double ey = m_y + y_k(0)*st + y_k(1)*ct - inv_s*x_k(1) - m_x*sx - m_y*cx;

        (*residuals)[k] = ex*ex + ey*ey;
    }
}
