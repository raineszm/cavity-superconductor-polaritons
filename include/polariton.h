#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <boost/math/tools/roots.hpp>
#include <cmath>
#include <complex>
#include <gsl/gsl_poly.h>

using boost::math::tools::toms748_solve;
using Eigen::Matrix3cd;
using Eigen::Vector3d;

#include "bs.h"
#include "cavity.h"
#include "coupling.h"

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

  Polariton(const BS& bs_, const Cavity& cav_, const Coupling& c_)
    : bs(bs_)
    , cav(cav_)
    , coupling(c_)
  {}

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
    auto par_factor = std::sqrt(2 / L);
    std::complex<double> c(0., par_factor * coupling.ImDA(omega));
    auto se = par_factor * par_factor * coupling.photon_se(omega, qx, qy);
    Matrix3cd act;
    act(0, 0) = bs.action(omega);
    act(0, 1) = c;
    act(1, 0) = -c;
    act.bottomRightCorner(2, 2) =
      cav.action(omega, qx, qy, coupling.state.sys.theta_v) + se;
    return act;
  }

  /** The eigenvalues of action()
   */
  Vector3d eigval(double omega, double qx, double qy) const
  {
    return action(omega, qx, qy).selfadjointView<Eigen::Upper>().eigenvalues();
  }

  /** Find a zero of action()
   */
  double find_mode(double qx, double qy) const
  {
    boost::uintmax_t max = 1e5;
    auto [a, b] = toms748_solve(
      [this, qx, qy](double omega) {
        return std::real(action(omega, qx, qy).determinant());
      },
      1e-3 * coupling.state.delta,
      1.99 * coupling.state.delta,
      [this, qx, qy](double a, double b) {
        double x = (a + b) / 2;
        return std::abs(action(x, qx, qy).determinant()) <
               1e-6 * coupling.state.delta;
      },
      max);
    return (a + b) / 2;
  }
};