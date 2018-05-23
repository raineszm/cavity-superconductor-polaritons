#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <array>
#include <boost/math/tools/roots.hpp>
#include <cmath>
#include <complex>
#include <gsl/gsl_poly.h>
#include <tuple>
#include <vector>

using boost::math::tools::newton_raphson_iterate;
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
    std::complex<double> c(0., std::sqrt(2 / L) * coupling.ImDA(omega));
    auto se = 2 / L * coupling.photon_se(omega, qx, qy);

    Matrix3cd act = Matrix3cd::Zero();
    act(0, 0) = bs.action(omega);
    act(0, 1) = c;
    act(1, 0) = -c;
    act.bottomRightCorner<2, 2>() =
      (cav.action(omega, qx, qy, coupling.state.sys.theta_v) + se)
        .cast<std::complex<double>>();
    return act;
  }

  Matrix3cd d_action(double omega, double qx, double qy) const
  {
    auto L = M_PI * C / cav.omega0;
    std::complex<double> c(0., std::sqrt(2 / L) * coupling.d_ImDA(omega));
    auto d_se = 2 / L * coupling.d_photon_se(omega, qx, qy);

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
    auto g = [this, qx, qy](double omega) { return det_and_d(omega, qx, qy); };

    std::array<double, 3> roots;
    double xl = 1e-3 * bs.root();
    double xu = 1.99 * coupling.state.delta;
    const double DX = 1e-3 * bs.root();

    boost::uintmax_t max = 1e5;
    roots[0] = newton_raphson_iterate(g, bs.root(), xl, xu, 32, max);

    std::vector<std::tuple<double, double>> intervals = {
      { xl, roots[0] - DX }, { roots[0] + DX, xu }
    };

    int count = 1;

    while (count < 3 && !intervals.empty()) {
      boost::uintmax_t max = 1e5;
      auto [xl, xu] = intervals.back();
      intervals.pop_back();
      auto root = newton_raphson_iterate(g, (xl + xu) / 2, xl, xu, 32, max);

      if (max != 1e5) {
        roots[count] = root;
        ++count;
        intervals.emplace_back(root + DX, xu);
        intervals.emplace_back(xl, root - DX);
      }
    }

    return roots;
  }
};