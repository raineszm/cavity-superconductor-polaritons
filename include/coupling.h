/**
 * @brief The coupling through the electronic system
 *
 * @file coupling.h
 * @author Zach Raines
 * @date 2018-06-12
 */
#pragma once
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <array>
#include <cmath>

#include "bs.h"
#include "cavity.h"
#include "state.h"

using Eigen::Matrix2d;

/**
 * @brief The coupling between the Superconductor and Cavity
 */
class Coupling
{
public:
  /**
   * @brief The Bardasis-Schrieffer mode
   *
   */
  const BS bs;
  /**
   * @brief The cavity photons
   *
   */
  const Cavity cav;

  /**
   * @brief Construct a new Coupling object
   *
   * @param bs_ #bs
   * @param cav_ #cav
   */
  explicit Coupling(const BS& bs_, const Cavity& cav_)
    : bs(bs_)
    , cav(cav_)
  {}

  /**
   * @cond
   */
  const State& state() const { return bs.state; }
  /**
   * @endcond
   */

  /**
   * @name Field Coupling
   *
   * Coupling between the field amplitudes \f$A,d\f$
   *
   * @{
   */

  /**
   * @brief The integrand for the imaginary part of the coupling
   *
   * @param l the BdG quasiparticle energy
   * @param theta the angle along the gap minumum
   * @param omega the frequency
   * @return double
   * @see ImDA()
   *
   * In evaluating the integral we pull out the constant prefactor
   * \f$-2 \frac{e}{c} \nu v_s i \Omega\Delta\f$
   * leaving the integrand
   * \f[
   * \frac{1}{\sqrt{\lambda^2 - \Delta^2}}
   * \frac{n_F(E^-(\lambda))-n_F(E^+(\lambda))}{(i\Omega_m)^2 -
   * (2\lambda)^2}f_d(\theta) \f]
   *
   * @note This result will appear in the retarded Green's function so we
   * analytically continue \f$i\Omega_m \to \omega + i0^+\f$.
   */
  double ImDA_int(double l, double theta, double omega) const
  {
    double fd = std::sqrt(2) * std::cos(2 * theta);
    double drift = state().sys.drift_theta(state().sys.kf(), theta);
    double Ep = drift + l;
    double Em = drift - l;

    return fd *
           (std::tanh(Ep / (2 * state().T)) - std::tanh(Em / (2 * state().T))) /
           ((omega * omega - 4 * l * l) * 2 *
            std::sqrt(l * l - state().delta * state().delta));
  }

  /**
   * @brief The integrand for the derivative of the imaginary part of the
   * coupling w.r.t \f$\omega\f$
   *
   * @param l the BdG quasiparticle energy
   * @param theta the angle along the gap minumum
   * @param omega the frequency
   * @return double
   * @see d_ImDA(), ImDA()
   */
  double d_ImDA_int(double l, double theta, double omega) const
  {
    double fd = std::sqrt(2) * std::cos(2 * theta);
    double drift = state().sys.drift_theta(state().sys.kf(), theta);
    double Ep = drift + l;
    double Em = drift - l;

    return -omega * fd *
           (std::tanh(Ep / (2 * state().T)) - std::tanh(Em / (2 * state().T))) /
           (std::pow(omega * omega - 4 * l * l, 2) *
            std::sqrt(l * l - state().delta * state().delta));
  }

  /**
   * @brief The coupling between the Bardasis Schrieffer mode and photons
   * appearing in the inverse Green's function
   *
   * @param omega the frequency
   * @return double
   * @see ImDA_int(), Polariton::inv_gf()
   *
   * This enters the action
   *  \f[ -i\sum_q g(i\Omega_m) \left(A^\parallel_q \
   * \bar{d}_{\perp,q} - A^{\parallel\ast}_q d_{\perp,q}\right)\f]
   *
   * with \f[ g_(i\Omega_m)
   * = \frac{e}{c} v_s i\Omega_m\Delta \sum_{\mathbf{k}}
   * \frac{n_F(E^-_\mathbf{k})-n_F(E^+_\mathbf{k}) }{(i\Omega_m)^2 -
   * (2\lambda_k)^2}\frac{f_d(\mathbf{k})}{\lambda_k} \f]
   * where here we use unrationalized (Gaussian) units
   *
   * We rewrite the integral in the \f$\xi\f$ approximation
   * \f[
   *  g(i\Omega_m) \approx 2 \frac{e}{c} \nu v_s i\Omega_m\Delta
   * \int_\Delta^\infty \frac{\lambda d\lambda}{\sqrt{\lambda^2 - \Delta^2}}
   * \int_0^{2\pi}\frac{d\theta}{2\pi}
   * \frac{n_F(E^-(\lambda))-n_F(E^+(\lambda))}{(i\Omega_m)^2 -
   * (2\lambda)^2}\frac{f_d(\theta)}{\lambda} \f]
   *
   * @note We take into account the angular factor arising from \f$v_s\cdot A\f$
   * when building the polariton matrix.
   *
   * @note This term will be analitically continued and appear in the inverse
   * inv_gf as \f$-g(\omega + i0^+)\f$
   *
   */
  double ImDA(double omega) const
  {
    return -2 * omega * state().delta * GPAR * state().sys.vs *
           state().sys.dos() *
           gsl_xi_integrate(
             [this, omega](double l, double theta) {
               return ImDA_int(l, theta, omega);
             },
             state().delta);
  }

  /**
   * @brief The derivative of the imaginary part of the coupling
   * \f$dg/d\omega\f$
   *
   * @param omega frequency
   * @return double
   * @see ImDA()
   */
  double d_ImDA(double omega) const
  {
    return -2 * state().delta * GPAR * state().sys.vs * state().sys.dos() *
           (gsl_xi_integrate(
              [this, omega](double l, double theta) {
                return ImDA_int(l, theta, omega);
              },
              state().delta) +
            omega * gsl_xi_integrate(
                      [this, omega](double l, double theta) {
                        return d_ImDA_int(l, theta, omega);
                      },
                      state().delta));
  }

  /** @}
   */

  /**
   * @name Photon Self-Energy
   *
   * @{
   */

  /**
   * @brief The photon self-energy due to renormalization by the s-wave state
   *
   * @param kx x component of electron momentum
   * @param ky y component of electron momentum
   * @param omega photon frequency
   * @param qx x component of photon momentum
   * @param qy y component of photon momentum
   * @param i left index
   * @param j right index
   * @param deriv is this the self enregy or its derivative?
   * @return double
   * @see State::, photon_se(), Polarition::inv_gf()
   *
   *
   * This term can be written in the form
   * \f[
   * \Pi^{ij}(\mathbf{q}, \Omega) =
   * \frac{e^2}{2c^2} \sum_\mathbf{k}
   * \sum_l \operatorname{tr}\left[\hat\pi_0(\mathbf{k} + \mathbf{q}/2,
   * \mathbf{k} - \mathbf{q}/2, \Omega) \hat\tau_l\right]T^{ij}_l(\mathbf{k},
   * \mathbf{q}) \f]
   *
   * where \f$\pi_0\f$ is given by @ref pi0().
   *
   * In this notation \f[
   * T^{ij}_0 =\ell_\mathbf{k,q}^2 v^i v^j + n_\mathbf{k,q}^2\, v_s^i v_s^j\\
   * T^{ij}_1 = p_\mathbf{k,q}^2 v^i v^j + m_\mathbf{k,q}^2\, v_s^i v_s^j\\
   * T^{ij}_2 = -i p_\mathbf{k,q}m_\mathbf{k,q}\left(v_s^i v^j + v^i
   * v_s^j\right)\\
   * T^{ij}_3 = \ell_\mathbf{k,q}n_\mathbf{k,q}\left(v_s^i v^j +
   * v^i v_s^j\right) \f]
   *
   * For covenience we take the the components \f$i=\{0,1\}\f$ to be along and
   * perpendicular to \f$\mathbf{v}_s\f$ respectively. Specifically, if
   * \f$\theta_s\f$ is the angle that \f$\mathbf{v}_s\f$ makes with the
   * \f$x\f$-axis, we define the new photon operators \f[ \begin{pmatrix}
   * A_\parallel\\
   * A_\perp
   * \end{pmatrix}
   * = \hat{R}(\theta_s)
   * \begin{pmatrix}
   * A_x\\
   * A_y
   * \end{pmatrix}
   * \f]
   * where \f$\hat{R}\f$ is a rotation matrix.
   * This induces the transformation \f$\mathbf{v} \to
   * R^{-1}(\theta_s)\mathbf{v}\f$ on the velocities. Concretely then \f[
   * \mathbf{v}_s &\to (v_s, 0)^T\\
   * \mathbf{v} &\to \frac{k}{m} (\cos(\theta - \theta_s),
   * \sin(\theta-\theta_s))^T \f]
   *
   */
  double photon_se_int(double kx,
                       double ky,
                       double omega,
                       double qx,
                       double qy,
                       int i,
                       int j,
                       bool deriv) const
  {
    const System& sys = state().sys;
    auto xp = sys.xi(kx + qx / 2, ky + qy / 2);
    auto xm = sys.xi(kx - qx / 2, ky - qy / 2);

    auto [k, theta_k] = gsl::rect_to_polar(kx, ky);
    auto vi = sys.v_comp(k, theta_k, i);
    auto vj = sys.v_comp(k, theta_k, j);
    auto vsi = sys.vs_comp(i);
    auto vsj = sys.vs_comp(j);

    auto T0 = state().l2(xp, xm) * vi * vj + state().n2(xp, xm) * vsi * vsj;
    auto T1 = state().p2(xp, xm) * vi * vj + state().m2(xp, xm) * vsi * vsj;
    auto iT2 = -state().mp(xp, xm) * (vsi * vj + vi * vsj);
    auto T3 = state().ln(xp, xm) * (vsi * vj + vi * vsj);

    auto lp = std::hypot(xp, state().delta);
    auto lm = std::hypot(xm, state().delta);

    auto dp = sys.drift(kx + qx / 2, ky + qy / 2);
    auto dm = sys.drift(kx - qx / 2, ky - qy / 2);

    auto pl = pi0_elems(dp, dm, lp, lm, omega, deriv);

    return T0 * pl[0] + T1 * pl[1] + iT2 * pl[2] + T3 * pl[3];
  }

  /**
   * @brief The polarization bubble entering into the photon self energy
   *
   * @param E1 quasiparticle energy
   * @param E2 quasiparticle energy
   * @param omega photon frequency
   * @param deriv is this the self-energy or its derivative w.r.t Freq
   * @return double
   * @see photon_se_int(), photon_se()
   *
   * \f[\pi_0(E_1, E_2, \omega) =
   * \frac{n_F(E_2) - n_F(E_1)}{i\Omega_m - E_1 + E_2}\f]
   *
   * For numerical convenience we make use of the relation
   * \f$n_f = \tfrac{1 - \tanh}{2}\f$ to rewrite this as
   * \f[\pi_0(E_1, E_2, \omega) =
   * \frac{1}{2}\frac{\tanh\frac{E_1}{2T}- \tanh\frac{E_2}{2T}}{i\Omega_m-
   * E_1 + E_2}\f]
   *
   * @note This is analytically continued to real frequency
   */
  double pi0(double E1, double E2, double omega, bool deriv) const
  {
    if (deriv) {
      return -0.5 *
             (std::tanh(E1 / (2 * state().T)) -
              std::tanh(E2 / (2 * state().T))) /
             std::pow(omega - E1 + E2, 2);

    } else {
      return 0.5 *
             (std::tanh(E1 / (2 * state().T)) -
              std::tanh(E2 / (2 * state().T))) /
             (omega - E1 + E2);
    }
  }

  /**
   * @brief The traces of pi0() with the pauli matrices.
   *
   * @param d1 doppler shift
   * @param d2 doppler shift
   * @param l1 quasiparticle energy
   * @param l2 quasi particle energy
   * @param omega photon frequency
   * @param deriv is this for the self-energy or its derivativer
   * @return std::array<double, 4>
   * @see pi0(), photon_se()
   *
   * Specifically we return
   * \f[\{\operatorname{tr}\pi_0\tau_0, \operatorname{tr}\pi_0\tau_1,
   * -i\operatorname{tr}\pi_0\tau_2, \operatorname{tr}\pi_0\tau_3\}\f]
   */
  std::array<double, 4> pi0_elems(double d1,
                                  double d2,
                                  double l1,
                                  double l2,
                                  double omega,
                                  bool deriv) const
  {
    auto p00 = pi0(d1 + l1, d2 + l2, omega, deriv);
    auto p01 = pi0(d1 + l1, d2 - l2, omega, deriv);
    auto p10 = pi0(d1 - l1, d2 + l2, omega, deriv);
    auto p11 = pi0(d1 - l1, d2 - l2, omega, deriv);

    return { { p00 + p11, p01 + p10, p01 - p10, p00 - p11 } };
  }

  //! @cond
  double _photon_se_or_deriv(double omega,
                             double qx,
                             double qy,
                             int i,
                             int j,
                             bool deriv) const
  {
    const auto& sys = state().sys;
    auto ret = gsl_angular_integrate(
      [this, omega, qx, qy, i, j, deriv](double k, double theta) {
        return photon_se_int(
          k * std::cos(theta), k * std::sin(theta), omega, qx, qy, i, j, deriv);
      },
      sys.kf() - 2 * state().delta / sys.vf(),
      sys.kf() + 2 * state().delta / sys.vf());

    // Diamagnetic?
    // if (i == j && !deriv) {
    //   ret += state().sys.n() / state().sys.m;
    // }

    return 0.5 * GPAR * GPAR * ret;
  }

  //!@endcond

  /**
   * @brief Evaluates the photon self-energy
   *
   * @param omega photon frequency
   * @param qx component of photon momentum
   * @param qy component of photon momentum
   * @param i left index of self-energy
   * @param j right index of self-energy
   * @return double
   * @see photon_se_int(), System::xi_k(), System::n()
   *
   * The self-energy is given by
   * \f[
   * \Pi^{ij}(\mathbf{q}, \Omega) =
   * \frac{e^2}{2c^2} \sum_\mathbf{k}
   * \sum_l \operatorname{tr}\left[\hat\pi_0(\mathbf{k} + \mathbf{q}/2,
   * \mathbf{k} - \mathbf{q}/2, \Omega) \hat\tau_l\right]T^{ij}_l(\mathbf{k},
   * \mathbf{q}) \f]
   * as described in photon_se_int().
   *
   * Numerically it is convenient to evaluate this integral in polar coordinates
   * via the change of variables \f$k=\sqrt{2m(\xi + \mu - \frac{1}{2}m
   * v_s^2)}\f$
   *
   * \f[ \Pi^{ij}(\mathbf{q}, \Omega) =
   *  \frac{e^2}{2c^2} \nu \int_{-\mu + \frac{1}{2}m
   * v_s^2}^\infty d\xi \int_0^{2\pi}\frac{d\theta}{2\pi} \sum_l
   * \operatorname{tr}\left[\hat\pi_0(\mathbf{k} + \mathbf{q}/2, \mathbf{k} -
   * \mathbf{q}/2, \Omega) \hat\tau_l\right]T^{ij}_l(\mathbf{k}, \mathbf{q}) \f]
   *
   * We then make use of the fact that \f$\mu\f$ is the largest energy scale in
   * the problem and take the lower limit of the \f$\xi\f$ integral to
   * \f$-\infty\f$ and symmetrize in \f$\xi\f$. Defining \f[ g(\xi, \theta,
   * \mathbf{q}) = \sum_l \operatorname{tr}\left[\hat\pi_0(\mathbf{k} +
   * \mathbf{q}/2, \mathbf{k} - \mathbf{q}/2, \Omega)
   * \hat\tau_l\right]T^{ij}_l(\mathbf{k}, \mathbf{q}) \f] i.e. photon_se_int(),
   * we can express the integral as
   *
   * \f[ \Pi^{ij}(\mathbf{q}, \Omega)\approx
   * \frac{e^2}{2c^2} \nu \int_0^\infty d\xi
   * \int_0^{2\pi}\frac{d\theta}{2\pi} \left[ g(\xi, \theta, \mathbf{q}) +
   * g(-\xi, \theta, \mathbf{q})\right] \f]
   */
  Matrix2d photon_se(double omega, double qx, double qy) const
  {
    return (Matrix2d() << _photon_se_or_deriv(omega, qx, qy, 0, 0, false),
            _photon_se_or_deriv(omega, qx, qy, 0, 1, false),
            _photon_se_or_deriv(omega, qx, qy, 1, 0, false),
            _photon_se_or_deriv(omega, qx, qy, 1, 1, false))
      .finished();
  }

  /**
   * @brief The derivative of the photon self energy w.r.t frequency
   *
   * @param omega photon frequency
   * @param qx coponenet of photon momentum
   * @param qy coponenet of photon momentum
   * @return Matrix2d
   *
   * \f[\frac{d\Pi(\omega, \mathbf{q})}{d\omega}\f]
   */
  Matrix2d d_photon_se(double omega, double qx, double qy) const
  {
    return (Matrix2d() << _photon_se_or_deriv(omega, qx, qy, 0, 0, true),
            _photon_se_or_deriv(omega, qx, qy, 0, 1, true),
            _photon_se_or_deriv(omega, qx, qy, 1, 0, true),
            _photon_se_or_deriv(omega, qx, qy, 1, 1, true))
      .finished();
  }

  //! @}

  /**
   * @name Mode Basis Coupling
   *
   * Coupling between the mode operators \f$a^\dagger,b^\dagger\f$
   * @{
   */

  /**
   * @brief The coupling between the Bardasis Schrieffer mode operator and
   * photon mode operator
   *
   * @param omega the frequency of the modes
   * @param qx x component of photon momentum
   * @param qy y component of photon momentum
   * @param i which photon polarization
   * @return double
   */
  double mode_coupling(double omega, double qx, double qy, int i) const
  {
    assert(i == 0 || i == 1);
    auto M = bs.M();
    auto V = cav.polarizations(qx, qy, state().sys.theta_s);
    auto vsdoteps = V(0, i);
    return 2 *
           std::sqrt(M_PI * C / (ALPHA * M * bs.root() * cav.omega(qx, qy))) *
           vsdoteps * ImDA(omega);
  }

  /**
   * @brief The derivative of the mode coupling \f$dg/d\omega\f$
   *
   * @param omega the frequency of the modes
   * @param qx x component of photon momentum
   * @param qy y component of photon momentum
   * @param i which photon polarization
   * @return double
   */
  double d_mode_coupling(double omega, double qx, double qy, int i) const
  {
    assert(i == 0 || i == 1);
    auto M = bs.M();
    auto V = cav.polarizations(qx, qy, state().sys.theta_s);
    auto vsdoteps = V(0, i);
    return 2 *
           std::sqrt(M_PI * C / (ALPHA * M * bs.root() * cav.omega(qx, qy))) *
           vsdoteps * d_ImDA(omega);
  }

  /**
   * @}
   */

  /**
   * @name Photon Mode Self Energy
   * @{
   */

  /**
   * @brief The self-energy as experienced by the photon mode operators
   *
   * @param omega photon frequency
   * @param qx component of momentum
   * @param qy component of momentum
   * @return Matrix2d
   *
   * The action is
   *
   * \f[
   *  S = \frac{1}{\beta}\sum_q \bar{a}_q \left[-i\omega_m + \omega_q +
   * \tilde{\Pi}(q)\right] a_q \f]
   *
   * In particular, (c.f. \ref photon_se(), \ref Polariton::inv_gf()), the
   * self energy we have previously defined enters as
   *
   * \f[
   *  \frac{-1}{\beta}\sum_q A_{-q} \Pi_q A_q
   * \f]
   *
   * that being the case the mode self-energy \f$\tilde{\Pi}\f$
   * is defined by
   *
   * \f[
   * \frac{2\pi c^2}{\omega_q} \epsilon(L/2)^* \Pi(q) \epsilon(L/2)
   * \f]
   * Our code implements the polarization vectors as \f$i\sqrt{\tfrac{2}{L}}V\f$
   * where $V$ is a real matrix describing the polarization vectors (c.f.
   * \ref Cavity::polarizations()) so the final expresion is
   *
   * \f[
   * \tilde{\Pi} = \frac{2\pi c^2}{\omega_q}\frac{2}{L} V^T \Pi(q) V
   * \f]
   */
  Matrix2d photon_se_mode(double omega, double qx, double qy) const
  {
    Matrix2d V = cav.polarizations(qx, qy, state().sys.theta_s);

    return 4 * M_PI * C / (ALPHA * cav.L() * cav.omega(qx, qy)) *
           V.transpose() *
           (Matrix2d() << _photon_se_or_deriv(omega, qx, qy, 0, 0, false),
            _photon_se_or_deriv(omega, qx, qy, 0, 1, false),
            _photon_se_or_deriv(omega, qx, qy, 1, 0, false),
            _photon_se_or_deriv(omega, qx, qy, 1, 1, false))
             .finished() *
           V;
  }

  /**
   * @brief The derivative of the mode self energy w.r.t frequency
   *
   * @param omega photon frequency
   * @param qx component of momentum
   * @param qy component of momentum
   * @return Matrix2d
   */
  Matrix2d d_photon_se_mode(double omega, double qx, double qy) const
  {
    Matrix2d V = cav.polarizations(qx, qy, state().sys.theta_s);
    return 4 * M_PI * C / (ALPHA * cav.L() * cav.omega(qx, qy)) *
           V.transpose() *
           (Matrix2d() << _photon_se_or_deriv(omega, qx, qy, 0, 0, true),
            _photon_se_or_deriv(omega, qx, qy, 0, 1, true),
            _photon_se_or_deriv(omega, qx, qy, 1, 0, true),
            _photon_se_or_deriv(omega, qx, qy, 1, 1, true))
             .finished() *
           V;
  }

  auto wf_renorm(double qx, double qy, double paraX) const
  {
    auto Z = Matrix2d::Identity() -
             paraX * paraX * d_photon_se_mode(cav.omega0, qx, qy);
    return Z.llt();
  }

  /**
   * @}
   */
};
