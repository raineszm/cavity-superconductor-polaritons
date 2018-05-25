#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <array>
#include <cmath>
#include <complex>
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
    auto se = 2 / L * big * big * coupling.photon_se(omega, qx, qy);

    return cav.action(omega, qx, qy, coupling.state.sys.theta_v) + se;
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

    Matrix3cd ret = Matrix3cd::Zero();
    ret(0, 0) = bs.d_action(omega);
    ret(0, 1) = c;
    ret(1, 0) = -c;
    ret.bottomRightCorner<2, 2>() =
      photon_sector(omega, qx, qy).cast<std::complex<double>>();
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
    auto detphot = [this, qx, qy](double omega) {
      return photon_sector(omega, qx, qy).determinant();
    };

    auto phot_perp = [this, qx, qy](double omega) {
      return std::real(action(omega, qx, qy)(2, 2));
    };

    auto other = [this, qx, qy](double omega) {
      auto act = action(omega, qx, qy);
      return std::real(act(0, 0) -
                       act(0, 1) * act(1, 0) * act(2, 2) /
                         act.bottomRightCorner<2, 2>().determinant());
    };
    std::array<double, 3> roots{ { NAN, NAN, NAN } };
    double xl = 0.5 * bs.root();
    double xu = 2 * bs.root();

    auto gsl_f = gsl_function_pp(phot_perp);
    auto solver = FSolver::create(gsl_root_fsolver_brent, gsl_f, xl, xu);

    auto perp_0 = solver.solve(1e-8 * bs.root(), 1e-6);

    if (gsl_root_test_residual(detphot(perp_0),
                               1e-8 * cav.omega0 * cav.omega0) == GSL_SUCCESS) {
      roots[0] = perp_0;
    }

    auto gsl_other = gsl_function_pp(other);
    solver.set(gsl_other, xl, xu);
    roots[1] = solver.solve(1e-8 * bs.root(), 1e-6);
    return roots;
  }
};