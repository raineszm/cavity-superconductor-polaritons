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
using Eigen::Matrix3d;
using Eigen::Vector2d;
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
   * @param c_ #coupling
   * @param paraX_ #paraX
   * @param dipoleX_ #dipoleX
   */
  Polariton(const Coupling& c_, double paraX_ = 1, double dipoleX_ = 1)
    : coupling(c_)
    , paraX(paraX_)
    , dipoleX(dipoleX_)
  {}

  virtual ~Polariton() = default;

  //! @cond
  const System& sys() const { return coupling.state().sys; }
  const State& state() const { return coupling.state(); }
  const Cavity& cav() const { return coupling.cav; }
  const BS& bs() const { return coupling.bs; }
  //! @endcond

  /**
   * @brief The inverse Green's function of the photon sector
   *
   * @param omega photon frequency
   * @param q photon momentum
   * @param theta_q angle of photon momentum w.r.t \f$v_s\f$
   * @return Matrix2d
   * @see Cavity::inv_gf(), Coupling::photon_se()
   *
   *
   * This includes the self energy term
   *
   * \f[
   * D^{-1}(i\Omega_m, \mathbf{q}) = D_0^{-1}((i\Omega_m)^2, \mathbf q) -
   * \Pi(i\Omega_m, \mathbf{q})\f]
   */
  Matrix2d photon_sector(double omega, double q, double theta_q) const
  {

    Matrix2d se = paraX * paraX * coupling.photon_se(omega, q, theta_q);
    Matrix2d cav_igf = cav().inv_gf(omega, q, theta_q);
    cav_igf -= se;
    return cav_igf;
  }

  /**
   * @brief Derivative of the photon sector w.r.t frequency
   * @param omega photon frequency
   * @param q component of photon momentum
   * @param theta_q component of photon momentum
   * @return Matrix2d
   * @see photon_sector()
   */
  Matrix2d d_photon_sector(double omega, double q, double theta_q) const
  {

    Matrix2d d_se = paraX * paraX * coupling.d_photon_se(omega, q, theta_q);
    Matrix2d d_cav_igf = cav().d_inv_gf(omega, q, theta_q);
    d_cav_igf -= d_se;
    return d_cav_igf;
  }

  /**
   * @brief The inverse inv_gf of the polariton
   *
   * @param omega frequency
   * @param q component of momentum
   * @param theta_q component of momentum
   * @return Matrix3cd
   * @see Coupling::ImDA(), Coupling::photon_se(), BS::inv_gf(),
   * Cavity::inv_gf(), photon_sector()
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
   * defines an inverse Greens' function which we analytically continue to real
   * frequency Here \f$D^{-1}\f$ is the bare cavity inverse inv_gf and
   * \f$G^{-1}\f$ is the BS inverse inv_gf. \f$A_\parallel\f$ and \f$A_\perp\f$
   * are the components of \f$\mathbf{A}\f$ parallel and perpendicular to the
   * supercurrent.
   */
  virtual Matrix3cd inv_gf(double omega, double q, double theta_q) const
  {
    std::complex<double> c(0., paraX * dipoleX * coupling.ImDA(omega));

    Matrix3cd act = Matrix3cd::Zero();
    act(0, 0) = bs().inv_gf(omega);
    act(0, 1) = c;
    act(1, 0) = -c;
    act.bottomRightCorner<2, 2>() =
      photon_sector(omega, q, theta_q).cast<std::complex<double>>();
    return act;
  }

  /**
   * @brief The determinant of the inverse Green's function of the fields
   *
   * @param omega frequnecy
   * @param q momentum
   * @param theta_q momentum
   * @return double
   * @see det_and_d()
   */
  double det_inv_gf(double omega, double q, double theta_q) const
  {
    return std::real(inv_gf(omega, q, theta_q).determinant());
  }

  /**
   * @brief The derivative of the inverse GF w.r.t frequency
   *
   * @param omega frequnecy
   * @param q momentum
   * @param theta_q momentum
   * @return Matrix3cd
   * @see inv_gf()
   */
  virtual Matrix3cd d_inv_gf(double omega, double q, double theta_q) const
  {
    std::complex<double> c(0., paraX * dipoleX * coupling.d_ImDA(omega));

    Matrix3cd act = Matrix3cd::Zero();
    act(0, 0) = bs().d_inv_gf(omega);
    act(0, 1) = c;
    act(1, 0) = -c;
    act.bottomRightCorner<2, 2>() =
      d_photon_sector(omega, q, theta_q).cast<std::complex<double>>();
    return act;
  }

  /**
   * @brief The determinant of the inverse GF and it's derivative
   *
   * @param omega frequency
   * @param q momentum
   * @param theta_q momentum
   * @return std::tuple<double, double> (det, deriv)
   * @see det_inv_gf(), d_det_inv_gf()
   */
  std::tuple<double, double> det_and_d(double omega,
                                       double q,
                                       double theta_q) const
  {

    auto A = inv_gf(omega, q, theta_q);
    double d_det =
      std::real((adjugate(A) * d_inv_gf(omega, q, theta_q)).trace());
    return std::tuple<double, double>(std::real(A.determinant()), d_det);
  }

  /**
   * @brief The derivative of the determinant of the inverse GF w.r.t frequency
   *
   * @param omega frequency
   * @param q momentum
   * @param theta_q momentum
   * @return double
   *
   * @see inv_gf(), d_inv_gf()
   */
  double d_det_inv_gf(double omega, double q, double theta_q) const
  {
    return std::get<1>(det_and_d(omega, q, theta_q));
  }

  /**
   * @brief The eigenvalues of inv_gf()
   *
   * @param omega frequency
   * @param q momentum
   * @param theta_q momentum
   * @return Vector3d
   */
  Vector3d eigval(double omega, double q, double theta_q) const
  {
    return inv_gf(omega, q, theta_q).eigenvalues().real();
  }

  /**
   * @brief Find the frequencies at which inv_gf() has a zero eigenvalue
   *
   * @param q momentum
   * @param theta_q momentum
   * @param ftol the tolerance in the derivative to be considered an extremum
   * @param double_root_tol the cutoff below which to declare something a double
   * root
   * @return std::array<double, 3> a sorted list of frequencies
   *
   * The frequencies returned correspond to the on-shell energies of the
   * associated modes.
   *
   * @note this returns a NaN if for those modes which are not found within the
   * chosen interval
   */
  std::array<double, 3> find_modes(double q,
                                   double theta_q,
                                   double ftol /* = 1e-10*/,
                                   double double_root_tol /* = 1e-17*/) const
  {

    // Define functions to be used in root finding
    auto det = [this, q, theta_q](double omega) {
      return det_inv_gf(omega, q, theta_q);
    };

    // Here we use the fact that our functions has two local extrema (since it
    // has 3 roots) We find the two extrem then look for the roots between each
    // and the outer boundary and the roots between them

    // Find the extrema
    auto gsl_det = make_gsl_function(det);

    std::array<double, 3> roots; // The return array
    roots.fill(std::numeric_limits<double>::quiet_NaN());

    const double xl = 0.6 * bs().root();
    const double xu = 1.9 * state().delta;
    // 1.99 * coupling.state.delta;
    auto extrema = _extrema(q, theta_q, xl, xu, ftol);

    FSolver fsolver(gsl_root_fsolver_brent);
    size_t count = 0;

    for (auto ext : extrema) {
      // Make sure we have an extremum
      if (std::isnan(ext)) {
        continue;
      }
      // Check for double root
      // We expect only one such double root
      if (std::fabs(det(ext)) < double_root_tol) {
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

  /**
   * @brief Find the local extrema of the determinant
   *
   * @param q momentum
   * @param theta_q momentum
   * @param xl lower frequency bound
   * @param xu upper frequency bound
   * @param ftol the tolerance in the derivative to be considered an extremum
   * @return std::array<double, 2> the positions of the extrema
   *
   * This is used in the rootfinding algorithm
   */
  std::array<double, 2> _extrema(double q,
                                 double theta_q,
                                 double xl,
                                 double xu,
                                 double ftol /*= 1e-10*/) const
  {

    auto d_det = [this, q, theta_q](double omega) {
      return d_det_inv_gf(omega, q, theta_q);
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

    extrema[0] = fdf_solver.solve(ftol, 100UL, EPS);

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

/**
 * @brief Find polaritons in the mode basis
 *
 */
class ModePolariton : public Polariton
{

public:
  using Polariton::Polariton;
  virtual ~ModePolariton() = default;
  /**
   * @brief The inverse inv_gf of the polariton
   *
   * @param omega frequency
   * @param q component of momentum
   * @param theta_q component of momentum
   * @return Matrix3cd
   * @see Coupling::mode_coupling(), Coupling::photon_se_mode(),
   * BS::hamiltonian(), Cavity::omega()
   *
   * The action
   * \f[
   * S_\text{pol} = -\frac{1}{\beta}\sum_q
   * \begin{pmatrix}
   * \bar{b}_q& \bar{a}_1(q)& \bar{a}_2(q)
   * \end{pmatrix}
   * \begin{pmatrix}
   * i\Omega - \Omega_\text{BS}&-g_\text{eff}(i\Omega_m, \mathbf{q})&0\\
   * -g_\text{eff}(i\Omega_m, \mathbf{q})& i\Omega_m - \omega_\mathbf{q} -
   * \Pi^{(0, 0)}(q)&
   * -\Pi^{(0, 1)}(q)\\ 0 & -\Pi^{(1, 0)}(q)& i \Omega_m - \omega_\mathbf{q}
   * \Pi^{(1,1)}(q) \end{pmatrix} \begin{pmatrix} b_q\\
   * a_1(q)\\ a_2(q) \end{pmatrix} \f]
   * defines an inverse Greens' function which we analytically continue to real
   * frequency
   */
  virtual Matrix3cd inv_gf(double omega,
                           double q,
                           double theta_q) const override
  {
    double c0 = paraX * dipoleX * coupling.mode_coupling(omega, q, theta_q, 0);
    double c1 = paraX * dipoleX * coupling.mode_coupling(omega, q, theta_q, 1);

    Matrix3d igf = Matrix3d::Zero();
    igf(0, 0) = omega - bs().hamiltonian();
    igf(0, 1) = -c0;
    igf(1, 0) = -c0;
    igf(0, 2) = -c1;
    igf(2, 0) = -c1;
    igf.bottomRightCorner<2, 2>() =
      (omega - cav().omega(q)) * Matrix2d::Identity() -
      paraX * paraX * coupling.photon_se_mode(omega, q, theta_q);
    return igf.cast<std::complex<double>>();
  }

  /**
   * @brief The derivative of the inverse GF w.r.t frequency
   *
   * @param omega frequnecy
   * @param q momentum
   * @param theta_q momentum
   * @return Matrix3cd
   * @see inv_gf()
   */
  virtual Matrix3cd d_inv_gf(double omega,
                             double q,
                             double theta_q) const override
  {
    double c0 =
      paraX * dipoleX * coupling.d_mode_coupling(omega, q, theta_q, 0);
    double c1 =
      paraX * dipoleX * coupling.d_mode_coupling(omega, q, theta_q, 1);
    Matrix3d igf = Matrix3d::Zero();
    igf(0, 0) = 1.;
    igf(0, 1) = -c0;
    igf(1, 0) = -c0;
    igf(0, 2) = -c1;
    igf(2, 0) = -c1;
    igf.bottomRightCorner<2, 2>() =
      Matrix2d::Identity() -
      paraX * paraX * coupling.d_photon_se_mode(omega, q, theta_q);
    return igf.cast<std::complex<double>>();
  }

  /**
   * @brief The hamiltonian of the coupled system
   *
   * @param q momentum
   * @param theta_q angle w.r.t \f$\mathbf{v}_s\f$
   * @return Matrix3d
   */
  Matrix3d hamiltonian(double q, double theta_q) const
  {
    Vector2d g{ coupling.mode_coupling(bs().root(), q, theta_q, 0),
                coupling.mode_coupling(bs().root(), q, theta_q, 1) };
    g *= paraX * dipoleX;

    auto LLT = coupling.wf_renorm(q, theta_q, paraX);
    auto L = LLT.matrixL();
    auto L_ = L.adjoint();
    L.solveInPlace(g);

    Matrix3d H = Matrix3d::Zero();
    H(0, 0) = bs().hamiltonian();
    H.block<1, 2>(0, 1) = g.adjoint();
    H.block<2, 1>(1, 0) = g;
    H.bottomRightCorner<2, 2>() =
      cav().omega(q) * Matrix2d::Identity() +
      paraX * paraX * coupling.photon_se_mode(cav().omega0, q, theta_q);
    H.bottomRightCorner<2, 2>() = L_.solve<Eigen::OnTheRight>(L.solve(
      cav().omega(q) * Matrix2d::Identity() +
      paraX * paraX * coupling.photon_se_mode(cav().omega0, q, theta_q)));
    return H;
  }
};
