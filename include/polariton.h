/**
 * @brief the polaritons
 *
 * @file polariton.h
 * @author Zach Raines
 * @date 2018-06-12
 */
#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <gsl/gsl_deriv.h>
#include <tuple>
#include <vector>

using Eigen::Matrix2d;
using Eigen::Matrix3cd;
using Eigen::Vector3d;

#include "bs.h"
#include "cavity.h"
#include "coupling.h"
#include "roots.h"
#include "utils.h"

/**
 * @brief The hybridized photon-Bardasis-Schrieffer object
 *
 */
class Polariton
{
public:
  //! The fermionic contribution/coupling
  const Coupling coupling;
  //! Enhancement of the paramagnetic coupling between electrons and the A field
  double paraX;
  //! Enhancment of the dipole coupling \f$g\f$ between the modes
  double dipoleX;

  /**
   * @brief Construct a new Polariton object
   *
   * @param c_
   * @param paraX_
   * @param dipoleX_
   */
  Polariton(const Coupling& c_, double paraX_ = 1, double dipoleX_ = 1)
    : coupling(c_)
    , paraX(paraX_)
    , dipoleX(dipoleX_)
  {}

  const System& sys() const { return coupling.state().sys; }
  const State& state() const { return coupling.state(); }
  const Cavity& cav() const { return coupling.cav; }
  const BS& bs() const { return coupling.bs; }

  Matrix2d photon_sector(double omega, double qx, double qy) const
  {

    Matrix2d se = paraX * paraX * coupling.photon_se(omega, qx, qy);
    Matrix2d cavaction = cav().gf(omega, qx, qy, sys().theta_v);
    cavaction += se;
    return cavaction;
  }

  Matrix2d d_photon_sector(double omega, double qx, double qy) const
  {

    Matrix2d d_se = paraX * paraX * coupling.d_photon_se(omega, qx, qy);
    Matrix2d d_cavaction = cav().d_gf(omega, qx, qy, sys().theta_v);
    d_cavaction += d_se;
    return d_cavaction;
  }

  /**
   * @brief The inverse GF of the polariton
   *
   * @param omega
   * @param qx
   * @param qy
   * @return Matrix3cd
   * @see Coupling::ImDA(), Coupling::photon_se(), BS::gf(), Cavity::gf()
   *
   * The action
   * \f[
   * S_\text{pol} = -\frac{1}{\beta}\sum_q
   * \begin{pmatrix}
   * d_\perp(-q)& A_\parallel(-q)& A_\perp(-q)
   * \end{pmatrix}
   * \begin{pmatrix}
   * S^{-1}(i\Omega_m)&ig(i\Omega_m)&0\\
   * -ig(i\Omega_m)& D_0^{-1}(q) - \Pi_{\parallel,\parallel}(q)&
   * -\Pi^{\parallel,\perp}(q)\\ 0 & -\Pi^{\perp,\parallel}(q)& D_0^{-1}(q) -
   * \Pi^{\perp,\perp}(q) \end{pmatrix} \begin{pmatrix} d_\perp(q)\\
   * A_\parallel(q)\\ A_\perp(q) \end{pmatrix} \f]
   * defines and inverse Greens' function which we analytically continue to real
   * frequency Here \f$D^{-1}\f$ is the bare cavity inverse GF and \f$G^{-1}\f$
   * is the BS inverse GF. \f$A_\parallel\f$ and \f$A_\perp\f$ are the
   * components of \f$\mathbf{A}\f$ parallel and perpendicular to the
   * supercurrent.
   */
  Matrix3cd gf(double omega, double qx, double qy) const
  {
    std::complex<double> c(0., paraX * dipoleX * coupling.ImDA(omega));

    Matrix3cd act = Matrix3cd::Zero();
    act(0, 0) = bs().gf(omega);
    act(0, 1) = c;
    act(1, 0) = -c;
    act.bottomRightCorner<2, 2>() =
      photon_sector(omega, qx, qy).cast<std::complex<double>>();
    return act;
  }

  double det_gf(double omega, double qx, double qy) const
  {
    return std::real(gf(omega, qx, qy).determinant());
  }

  /** The derivative of the inverse GF
   */
  Matrix3cd d_gf(double omega, double qx, double qy) const
  {
    std::complex<double> c(0., paraX * dipoleX * coupling.d_ImDA(omega));

    Matrix3cd act = Matrix3cd::Zero();
    act(0, 0) = bs().d_gf(omega);
    act(0, 1) = c;
    act(1, 0) = -c;
    act.bottomRightCorner<2, 2>() =
      d_photon_sector(omega, qx, qy).cast<std::complex<double>>();
    return act;
  }

  std::tuple<double, double> det_and_d(double omega, double qx, double qy) const
  {

    auto A = gf(omega, qx, qy);
    double d_det = std::real((adjugate(A) * d_gf(omega, qx, qy)).trace());
    return std::tuple<double, double>(std::real(A.determinant()), d_det);
  }

  double d_det_gf(double omega, double qx, double qy) const
  {
    return std::get<1>(det_and_d(omega, qx, qy));
  }

  /** The eigenvalues of gf()
   */
  Vector3d eigval(double omega, double qx, double qy) const
  {
    return gf(omega, qx, qy).eigenvalues().real();
  }

  /** Find zeroes of gf()
   */
  std::array<double, 3> find_modes(double qx, double qy) const
  {

    // Define functions to be used in root finding
    auto det = [this, qx, qy](double omega) { return det_gf(omega, qx, qy); };

    // Here we use the fact that our functions has two local extrema (since it
    // has 3 roots) We find the two extrem then look for the roots between each
    // and the outer boundary and the roots between them

    // Find the extrema
    auto gsl_det = make_gsl_function(det);

    std::array<double, 3> roots; // The return array
    roots.fill(std::numeric_limits<double>::quiet_NaN());

    const double xl = 0.6 * bs().root();
    const double xu = 1.5 * bs().root();
    // 1.99 * coupling.state.delta;
    auto extrema = _extrema(qx, qy, xl, xu);

    FSolver fsolver(gsl_root_fsolver_brent);
    size_t count = 0;

    for (auto ext : extrema) {
      // Make sure we have an extremum
      if (std::isnan(ext)) {
        continue;
      }
      // Check for double root
      // We expect only one such double root
      if (std::fabs(det(ext)) < 1e-17) {
        roots[count++] = ext;
        roots[count++] = ext;
      } else {
        try {
          // Find root beyond first ext
          if (det(xl) * det(ext) > 0) {
            fsolver.set(gsl_det, ext, xu);
          } else {
            fsolver.set(gsl_det, xl, ext);
          }

          roots[count++] = fsolver.solve(1e-8 * bs().root(), 1e-8);
        } catch (const gsl::GSLException&) {
        }
      }
    }

    if (std::isnan(extrema[0])) {
      extrema[0] = xl;
    }
    if (std::isnan(extrema[1])) {
      extrema[1] = xu;
    }

    try {
      // Find root between extrema
      fsolver.set(gsl_det, extrema[0], extrema[1]);
      roots[count++] = fsolver.solve(1e-8 * bs().root(), 1e-8);
    } catch (const gsl::GSLException&) {
    }

    // Check if we've had any double roots
    std::sort(roots.begin(), roots.end());

    return roots;
  }

  std::array<double, 2> _extrema(double qx,
                                 double qy,
                                 double xl,
                                 double xu) const
  {

    auto d_det = [this, qx, qy](double omega) {
      return d_det_gf(omega, qx, qy);
    };

    auto gsl_d_det = make_gsl_function(d_det);

    const double EPS = 1e-6 * bs().root();

    auto d_det_fdf = [EPS, gsl_d_det](double omega) {
      auto d_det_val = GSL_FN_EVAL(&gsl_d_det, omega);
      auto d_d_det = deriv_gsl(gsl_d_det, omega, EPS);
      return std::tuple<double, double>(d_det_val, d_d_det);
    };

    std::array<double, 2> extrema;
    extrema.fill(std::numeric_limits<double>::quiet_NaN());

    FDFSolver fdf_solver(gsl_root_fdfsolver_newton);
    auto gsl_d_det_fdf = make_gsl_function_fdf(d_det, d_det_fdf);
    fdf_solver.set(gsl_d_det_fdf, bs().root());

    extrema[0] = fdf_solver.solve(1e-10, 100UL, EPS);

    FSolver fsolver(gsl_root_fsolver_brent);

    try {
      if (deriv_gsl(gsl_d_det, extrema[0], EPS) * d_det(xu) > 0) {
        fsolver.set(gsl_d_det, xl, extrema[0] - EPS);
      } else {
        fsolver.set(gsl_d_det, extrema[0] + EPS, xu);
      }

      // Find other extremum
      extrema[1] = fsolver.solve(1e-8 * bs().root(), 1e-8);
    } catch (const gsl::GSLException&) {
      extrema[1] = std::numeric_limits<double>::quiet_NaN();
    }

    std::sort(extrema.begin(), extrema.end());
    return extrema;
  }
};
