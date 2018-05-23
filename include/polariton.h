#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <array>
#include <cmath>
#include <complex>
#include <gsl/gsl_errno.h>
#include <tuple>
#include <vector>

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
    auto se = 2 / L * big * big * coupling.photon_se(omega, qx, qy);

    Matrix3cd act = Matrix3cd::Zero();
    act(0, 0) = bs.action(omega);
    act(0, 1) = c;
    act(1, 0) = -c;
    act.bottomRightCorner<2, 2>() =
      (cav.action(omega, qx, qy, coupling.state.sys.theta_v) + se)
        .cast<std::complex<double>>();
    return act;
  }

  /** The derivative of the inverse GF
   */
  Matrix3cd d_action(double omega, double qx, double qy) const
  {
    auto L = M_PI * C / cav.omega0;
    std::complex<double> c(0., big * std::sqrt(2 / L) * coupling.d_ImDA(omega));
    auto d_se = 2 / L * big * big * coupling.d_photon_se(omega, qx, qy);

    Matrix3cd ret = Matrix3cd::Zero();
    ret(0, 0) = bs.d_action(omega);
    ret(0, 1) = c;
    ret(1, 0) = -c;
    ret.bottomRightCorner<2, 2>() =
      (cav.d_action(omega, qx, qy, coupling.state.sys.theta_v) + d_se)
        .cast<std::complex<double>>();
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
    return action(omega, qx, qy).selfadjointView<Eigen::Upper>().eigenvalues();
  }

  /** Find zeroes of action()
   */
  std::array<double, 3> find_modes(double qx, double qy) const
  {
    auto f = [this, qx, qy](double omega) {
      return std::real(action(omega, qx, qy).determinant());
    };
    auto g = [this, qx, qy](double omega) { return det_and_d(omega, qx, qy); };

    std::array<double, 3> roots;
    double xl = 0.5 * bs.root();
    double xu = 2 * bs.root();
    const double DX = 1e-3 * bs.root();

    auto solver = FDFSolver(gsl_root_fdfsolver_steffenson);

    auto gsl_fdf = gsl_function_fdf_pp(f, g);
    auto last = bs.root();
    solver.set(gsl_fdf, last);

    const int NITER = 100;

    for (int i = 0; i < NITER; i++) {
      solver.step();
      if (gsl_root_test_delta(solver.root(), last, 1e-3 * bs.root(), 1e-3) ==
          GSL_SUCCESS) {
        roots[0] = solver.root();
        break;
      }

      last = solver.step();
    }

    // last = 0.5 * bs.root();
    // solver.set(gsl_fdf, last);
    // for (int i = 0; i < NITER; i++) {
    //   solver.step();
    //   if (gsl_root_test_delta(solver.root(), last, 1e-3 * bs.root(), 1e-3) ==
    //       GSL_SUCCESS) {
    //     roots[1] = solver.root();
    //     break;
    //   }

    //   last = solver.step();
    // }

    last = 1.5 * bs.root();
    solver.set(gsl_fdf, last);
    for (int i = 0; i < NITER; i++) {
      solver.step();
      if (gsl_root_test_delta(solver.root(), last, 1e-3 * bs.root(), 1e-3) ==
          GSL_SUCCESS) {
        roots[2] = solver.root();
        break;
      }

      last = solver.step();
    }

    return roots;
  }
};