#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <tuple>
#include <vector>

using Eigen::Matrix2d;
using Eigen::Matrix3cd;
using Eigen::Vector3d;

#include "bs.h"
#include "cavity.h"
#include "coupling.h"
#include "roots.h"

/** The hybridized photon-Bardasis-Schrieffer object
 */
class Polariton
{
public:
  //! The Bardasis Schrieffer mode
  const BS bs;
  //! The photonic sector
  const Cavity cav;
  //! The fermionic contribution/coupling
  const Coupling coupling;
  const double big;

  Polariton(const BS& bs_,
            const Cavity& cav_,
            const Coupling& c_,
            double BIGUP = 1)
    : bs(bs_)
    , cav(cav_)
    , coupling(c_)
    , big(BIGUP)
  {
    if (bs.state != coupling.state) {
      throw std::invalid_argument(
        "bs and c must be constructed with the same state");
    }
  }

  Matrix2d photon_sector(double omega, double qx, double qy) const
  {

    auto L = M_PI * C / cav.omega0;
    Matrix2d se = 2 / L * big * big * coupling.photon_se(omega, qx, qy);
    Matrix2d cavaction = cav.action(omega, qx, qy, coupling.state.sys.theta_v);
    cavaction += se;
    return cavaction;
  }

  /** The inverse GF of the polariton
   * \f[
   * S_\text{pol} = \sum_q
   * \begin{pmatrix}
   * d_\perp(-q)& A_\parallel(-q)& A_\perp(-q)
   * \end{pmatrix}
   * \begin{pmatrix}
   * S^{-1}(\Omega)&ig(\Omega)&0\\
   * -ig(\Omega)& D_0^{-1}(q) - \Pi_{\parallel,\parallel}(q)&
   * -\Pi^{\parallel,\perp}(q)\\ 0 & -\Pi^{\perp,\parallel}(q)& D_0^{-1}(q) -
   * \Pi^{\perp,\perp}(q) \end{pmatrix} \begin{pmatrix} d_\perp(q)\\
   * A_\parallel(q)\\ A_\perp(q) \end{pmatrix} \f] This method calculates the
   * central matrix of the above. Here \f$D^{-1}\f$ is the bare cavity inverse
   * GF and \f$S^{-1}\f$ is the BS inverse GF.
   * \f$A_\parallel\f$ and \f$A_\perp\f$ are the components of \f$\mathbf{A}\f$
   * parallel and perpendicular to the supercurrent.
   *
   * \sa Coupling::ImDA(), Coupling::photon_se(), BS::action(),
   * Cavity::action();
   */
  Matrix3cd action(double omega, double qx, double qy) const
  {
    auto L = M_PI * C / cav.omega0;
    std::complex<double> c(0., big * std::sqrt(2 / L) * coupling.ImDA(omega));

    Matrix3cd act = Matrix3cd::Zero();
    act(0, 0) = bs.action(omega);
    act(0, 1) = c;
    act(1, 0) = -c;
    act.bottomRightCorner<2, 2>() =
      photon_sector(omega, qx, qy).cast<std::complex<double>>();
    return act;
  }

  /** The derivative of the inverse GF
   */
  Matrix3cd d_action(double omega, double qx, double qy) const
  {
    auto L = M_PI * C / cav.omega0;
    std::complex<double> c(0., big * std::sqrt(2 / L) * coupling.d_ImDA(omega));

    auto d_se = 2 / L * big * big * coupling.d_photon_se(omega, qx, qy);
    auto d_cavaction = cav.d_action(omega, qx, qy, coupling.state.sys.theta_v);
    d_cavaction += d_se;

    Matrix3cd ret = Matrix3cd::Zero();
    ret(0, 0) = bs.d_action(omega);
    ret(0, 1) = c;
    ret(1, 0) = -c;
    ret.bottomRightCorner<2, 2>() = d_cavaction.cast<std::complex<double>>();
    return ret;
  }

  std::tuple<double, double> det_and_d(double omega, double qx, double qy) const
  {

    auto A = action(omega, qx, qy);
    double d_det = std::real((adjugate(A) * d_action(omega, qx, qy)).trace());
    return std::tuple<double, double>(std::real(A.determinant()), d_det);
  }

  /** The eigenvalues of action()
   */
  Vector3d eigval(double omega, double qx, double qy) const
  {
    return action(omega, qx, qy).eigenvalues().real();
  }

  /** Find zeroes of action()
   */
  std::array<double, 3> find_modes(double qx, double qy) const
  {
    const double EPS = 1e-6 * bs.root();

    // Define functions to be used in root finding
    auto det = [this, qx, qy](double omega) {
      return std::real(action(omega, qx, qy).determinant());
    };

    auto d_det = [this, qx, qy](double omega) {
      return std::get<1>(det_and_d(omega, qx, qy));
    };

    gsl_function_pp<decltype(d_det)> gsl_d_det(d_det);

    auto d_det_fdf = [this, gsl_d_det](double omega) {
      auto d_det_val = GSL_FN_EVAL(&gsl_d_det, omega);
      double d_d_det, err;
      gsl_deriv_central(&gsl_d_det, omega, 1e-3 * bs.root(), &d_d_det, &err);
      return std::tuple<double, double>(d_det_val, d_d_det);
    };

    std::array<double, 3> roots; // The return array
    roots.fill(std::numeric_limits<double>::quiet_NaN());

    FDFSolver fdf_solver(gsl_root_fdfsolver_steffenson);

    gsl_function_fdf_pp<decltype(d_det), decltype(d_det_fdf)> gsl_d_det_fdf(
      d_det, d_det_fdf);

    fdf_solver.set(gsl_d_det_fdf, bs.root());

    // Here we use the fact that our functions has two local extrema (since it
    // has 3 roots) We find the two extrem then look for the roots between each
    // and the outer boundary and the roots between them

    // Find the first extremum
    double extremum = fdf_solver.solve(1e-18);

    double deriv2, err;
    const double xl = 0.5 * bs.root();
    const double xu = 1.99 * bs.root();
    FSolver fsolver(gsl_root_fsolver_brent);
    gsl_deriv_central(&gsl_d_det, extremum, EPS, &deriv2, &err);

    if (deriv2 * d_det(xu) > 0) {
      fsolver.set(gsl_d_det, xl, extremum - EPS);
    } else {
      fsolver.set(gsl_d_det, extremum + EPS, xu);
    }

    // Find other extremum
    auto extremum2 = fsolver.solve(1e-8 * bs.root(), 1e-8);

    size_t count = 0;
    gsl_function_pp<decltype(det)> gsl_det(det);

    // Check for double root
    // We expect only one such double root
    if (gsl_root_test_residual(det(extremum), 1e-20) == GSL_SUCCESS) {
      roots[count++] = extremum;
      roots[count++] = extremum;
    } else {
      // Find root beyond first extremum
      if (det(xl) * det(extremum) > 0) {
        fsolver.set(gsl_det, extremum, xu);
      } else {
        fsolver.set(gsl_det, xl, extremum);
      }

      roots[count++] = fsolver.solve(1e-8 * bs.root(), 1e-8);
    }

    // check again for double root
    if (gsl_root_test_residual(det(extremum2), 1e-20) == GSL_SUCCESS) {
      roots[count++] = extremum2;
      roots[count++] = extremum2;
    } else {
      // If not find root beyond second extremum
      if (det(xl) * det(extremum2) > 0) {
        fsolver.set(gsl_det, extremum2, xu);
      } else {
        fsolver.set(gsl_det, xl, extremum2);
      }
      roots[count++] = fsolver.solve(1e-8 * bs.root(), 1e-8);
    }

    // Find root between extrema
    if (extremum > extremum2) {
      fsolver.set(gsl_det, extremum2, extremum);
    } else {
      fsolver.set(gsl_det, extremum, extremum2);
    }
    roots[count++] = fsolver.solve(1e-8 * bs.root(), 1e-8);

    // Check if we've had any double roots
    std::sort(roots.begin(), roots.end());

    return roots;
  }
};
