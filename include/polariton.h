#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <boost/math/tools/roots.hpp>
#include <cmath>
#include <complex>

using boost::math::tools::toms748_solve;
using Eigen::Matrix3cd;
using Eigen::Vector3cd;

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
   * d_\perp(-q)& A^x(-q)& A^y(-q)
   * \end{pmatrix}
   * \begin{pmatrix}
   * S^{-1}(\Omega)&ig^x(\Omega)&ig^y(\Omega)\\
   * -ig^x(\Omega)& D_0^{-1}(q) - \Pi^{00}(q)& -\Pi^{01}(q)\\
   * -ig^y(\Omega)& -\Pi^{10}(q)& D_0^{-1}(q) - \Pi^{11}(q)
   * \end{pmatrix}
   * \begin{pmatrix}
   * d_\perp(q)\\ A^x(q)\\ A^y(q)
   * \end{pmatrix}
   * \f]
   * This method calculates the central matrix of the above.
   * Here \f$D^{-1}\f$ is the bare cavity inverse GF and \f$S^{-1}\f$ is the BS
   * inverse GF.
   *
   * \sa Coupling::ImDA(), Coupling::photon_se(), BS::action(),
   * Cavity::action();
   */
  Matrix3cd action(double omega, double qx, double qy) const
  {
    Matrix3cd mat;
    std::complex<double> c(0., coupling.ImDA(omega));
    auto cs = std::cos(coupling.state.sys.theta_v);
    auto sn = std::sin(coupling.state.sys.theta_v);
    auto se00 = coupling.photon_se(omega, qx, qy, 0, 0);
    auto se01 = coupling.photon_se(omega, qx, qy, 0, 1);
    auto se10 = coupling.photon_se(omega, qx, qy, 1, 0);
    auto se11 = coupling.photon_se(omega, qx, qy, 1, 1);
    mat << bs.action(omega), c * cs, c * sn, -c * cs,
      cav.action(omega, qx, qy) + se00, se01, -c * sn, se10,
      cav.action(omega, qx, qy) + se11;
    return mat;
  }

  /** The eigenvalues of action()
   */
  Vector3cd eigen(double omega, double qx, double qy) const
  {
    return action(omega, qx, qy).eigenvalues();
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