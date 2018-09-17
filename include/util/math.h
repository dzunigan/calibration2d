
// STL
#include <cmath>
#include <limits>

// Eigen
#include <Eigen/Core>
#include <Eigen/Dense>

template<typename Scalar>
inline bool is_negligible(Scalar x) {
    return std::abs(x) <= std::numeric_limits<Scalar>::epsilon();
}

/* Solve for real roots of the standard quadratic equation,
 * returning the number of real roots found.
 */
/* solve_quadratic.c - finds the real roots of a x^2 + b x + c = 0 */
// Adapted from gsl_poly.h (https://github.com/ampl/gsl/blob/master/poly/gsl_poly.h)
int solveQuadratic(double a, double b, double c, double &x0, double &x1) {
    if (is_negligible(a)) { /* Handle linear case */
        if (is_negligible(b))
            return 0;
        else {
          x0 = -c / b;
          return 1;
        }
    }

    double disc = b * b - 4 * a * c;
    if (disc > 0.0) {
        if (is_negligible(b)) {
            double r = std::sqrt(-c / a);
            x0 = -r;
            x1 =  r;
        } else {
            // What Every Computer Scientist Should Know About Floating-Point Arithmetic, Goldberg, D. in Computing Surveys, 1991 (Sec. 1.4 - Cancellation)
            // www.itu.dk/~sestoft/bachelor/IEEE754_article.pdf
            // Stable quadratic roots according to BKP Horn.
            // http://people.csail.mit.edu/bkph/articles/Quadratics.pdf
            double sgnb = (b > 0 ? 1 : -1);
            double temp = -0.5 * (b + sgnb * std::sqrt(disc));
            double r1 = temp / a ;
            double r2 = c / temp ;

            if (r1 < r2) {
                x0 = r1;
                x1 = r2;
            } else {
                x0 = r2;
                x1 = r1;
            }
        }

        return 2;
    } else if (is_negligible(disc)) {
        x0 = -0.5 * b / a;
        x1 = -0.5 * b / a;

        return 2;
    } else
        return 0;
}


// Solves the right nullspace from QR decomposition,
// returning the size of the kernel
template<typename Derived, typename Scalar>
int solveNullspace(const Eigen::MatrixBase<Derived> &A, Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &k) {
    Eigen::ColPivHouseholderQR<Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>> qr(A.transpose());
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> Q = qr.householderQ();

    int n = qr.dimensionOfKernel();
    k.resize(Q.rows(), n);

    k = Q.block(0, Q.cols() - n, Q.rows(), n);
    return qr.dimensionOfKernel();
}
