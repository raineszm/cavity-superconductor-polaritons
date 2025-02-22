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
    double fd = std::sqrt(2); // * std::cos(2 * (theta + state().sys.theta_s));
    double doppler = state().sys.doppler(state().sys.kf(), theta);
    double Ep = doppler + l;
    double Em = doppler - l;

    return fd * diff_of_tanh(Ep, Em, state().T) /
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
    double fd = std::sqrt(2); //* std::cos(2 * (theta + state().sys.theta_s));
    double doppler = state().sys.doppler(state().sys.kf(), theta);
    double Ep = doppler + l;
    double Em = doppler - l;

    return -omega * fd *
           (std::tanh(Ep / (2 * state().T)) - std::tanh(Em / (2 * state().T))) /
           (std::pow(omega * omega - 4 * l * l, 2) *
            std::sqrt(l * l - state().delta * state().delta));
  }

  /**
   * @brief The integral appearing in the definition of coupling
   *
   * @param omega frequency
   * @return double
   * @see ImDA(), mode_coupling()
   *
   *
   * \f[
   * 2\nu \int_\Delta^\infty d\lambda
   * \int_0^{2\pi}\frac{d\theta}{2\pi}\frac{1}{\sqrt{\lambda^2 - \Delta^2}}
   * \frac{n_F(E^-(\lambda))-n_F(E^+(\lambda))}{(i\Omega_m)^2 -
   * (2\lambda)^2}f_d(\theta) \f]
   */
  double coupling_I(double omega) const
  {
    return state().sys.dos() * gsl_xi_cos_integrate(
                                 [this, omega](double l, double theta) {
                                   return ImDA_int(l, theta, omega) -
                                          ImDA_int(l, M_PI_2 - theta, omega);
                                 },
                                 state().delta,
                                 state().sys.theta_s);
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
    return -omega * state().delta * GPAR * state().sys.vs * coupling_I(omega);
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
           (gsl_xi_cos_integrate(
              [this, omega](double l, double theta) {
                return ImDA_int(l, theta, omega);
              },
              state().delta,
              state().sys.theta_s) +
            omega * gsl_xi_cos_integrate(
                      [this, omega](double l, double theta) {
                        return d_ImDA_int(l, theta, omega);
                      },
                      state().delta,
                      state().sys.theta_s));
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
   * @param k electron momentum
   * @param theta_k electron angle
   * @param omega photon frequency
   * @param q photon momentum
   * @param theta_q photon angle
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
   * T^{ij}_1 = -p_\mathbf{k,q}^2 v^i v^j + m_\mathbf{k,q}^2\, v_s^i v_s^j\\
   * iT^{ij}_2 = p_\mathbf{k,q}m_\mathbf{k,q}\left(v^i v_s^j - v_s^i
   * v^j\right)\\
   * T^{ij}_3 = \ell_\mathbf{k,q}n_\mathbf{k,q}\left(v_s^i v^j +
   * v^i v_s^j\right) \f]
   *
   * @notes We take the \f$x\f$ axis to be defined by the supercurrent
   */
  double photon_se_int(double k,
                       double theta_k,
                       double omega,
                       double q,
                       double theta_q,
                       int i,
                       int j,
                       bool deriv) const
  {
    const System& sys = state().sys;
    auto theta_qk = theta_q - theta_k;
    auto kp = std::sqrt(k * k + q * q / 4 + k * q * std::cos(theta_qk));
    auto km = std::sqrt(k * k + q * q / 4 - k * q * std::cos(theta_qk));

    auto xp = sys.xi_k(kp);
    auto xm = sys.xi_k(km);

    auto vi = sys.v_comp(k, theta_k, i);
    auto vj = sys.v_comp(k, theta_k, j);
    auto vsi = sys.vs_comp(i);
    auto vsj = sys.vs_comp(j);

    auto T0 = state().l2(xp, xm) * vi * vj + state().n2(xp, xm) * vsi * vsj;
    auto T1 = -state().p2(xp, xm) * vi * vj + state().m2(xp, xm) * vsi * vsj;
    auto iT2 = state().mp(xp, xm) * (vi * vsj - vsi * vj);
    auto T3 = state().ln(xp, xm) * (vsi * vj + vi * vsj);

    auto lp = std::hypot(xp, state().delta);
    auto lm = std::hypot(xm, state().delta);

    auto theta_p = theta_k + std::atan2(q / 2 * std::sin(theta_qk),
                                        k + q / 2 * std::cos(theta_qk));
    auto theta_m = theta_k + std::atan2(-q / 2 * std::sin(theta_qk),
                                        k - q / 2 * std::cos(theta_qk));
    auto dp = sys.doppler(kp, theta_p);
    auto dm = sys.doppler(km, theta_m);

    auto pl = pi0_elems(dp, dm, lp, lm, omega, deriv);

    return T0 * pl[0] + T1 * pl[1] + iT2 * pl[2] + T3 * pl[3];
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
    const State& state_ = state();
    auto p00 = state_.pi0(d1 + l1, d2 + l2, omega, deriv);
    auto p01 = state_.pi0(d1 + l1, d2 - l2, omega, deriv);
    auto p10 = state_.pi0(d1 - l1, d2 + l2, omega, deriv);
    auto p11 = state_.pi0(d1 - l1, d2 - l2, omega, deriv);

    return { { p00 + p11, p01 + p10, p01 - p10, p00 - p11 } };
  }

  //! @cond
  double _photon_se_or_deriv(double omega,
                             double q,
                             double theta_q,
                             int i,
                             int j,
                             bool deriv) const
  {
    const auto& sys = state().sys;
    auto ret = gsl_angular_integrate(
      [this, omega, q, theta_q, i, j, deriv](double k, double theta) {
        return photon_se_int(k, theta, omega, q, theta_q, i, j, deriv);
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
   * @param q photon momentum
   * @param theta_q angle of photon momentum w.r.t \f$v_s\f$
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
   */
  Matrix2d photon_se(double omega, double q, double theta_q) const
  {
    return (Matrix2d() << _photon_se_or_deriv(omega, q, theta_q, 0, 0, false),
            _photon_se_or_deriv(omega, q, theta_q, 0, 1, false),
            _photon_se_or_deriv(omega, q, theta_q, 1, 0, false),
            _photon_se_or_deriv(omega, q, theta_q, 1, 1, false))
      .finished();
  }

  /**
   * @brief The derivative of the photon self energy w.r.t frequency
   *
   * @param omega photon frequency
   * @param q photon momentum
   * @param theta_q angle of photon momentum w.r.t \f$v_s\f$
   * @return Matrix2d
   *
   * \f[\frac{d\Pi(\omega, \mathbf{q})}{d\omega}\f]
   */
  Matrix2d d_photon_se(double omega, double q, double theta_q) const
  {
    return (Matrix2d() << _photon_se_or_deriv(omega, q, theta_q, 0, 0, true),
            _photon_se_or_deriv(omega, q, theta_q, 0, 1, true),
            _photon_se_or_deriv(omega, q, theta_q, 1, 0, true),
            _photon_se_or_deriv(omega, q, theta_q, 1, 1, true))
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
   * @param q photon momentum
   * @param theta_q angle of photon momentum w.r.t \f$v_s\f$
   * @param i which photon polarization
   * @return double
   */
  double mode_coupling(double omega, double q, double theta_q, int i) const
  {
    assert(i == 0 || i == 1);
    auto M = bs.M();
    auto V = cav.polarizations(q, theta_q);
    auto ivsdoteps = state().sys.vs * std::sqrt(2 / cav.L()) * V(0, i);
    auto I = coupling_I(omega);
    return -ivsdoteps * state().delta *
           std::sqrt((4 * M_PI * ALPHA * C * bs.root()) / (M * cav.omega(q))) *
           I;
  }

  /**
   * @brief The derivative of the mode coupling \f$dg/d\omega\f$
   *
   * @param omega the frequency of the modes
   * @param q photon momentum
   * @param theta_q angle of photon momentum w.r.t \f$v_s\f$
   * @param i which photon polarization
   * @return double
   */
  double d_mode_coupling(double omega, double q, double theta_q, int i) const
  {
    assert(i == 0 || i == 1);
    auto M = bs.M();
    auto V = cav.polarizations(q, theta_q);
    auto vsdoteps = V(0, i);
    return 2 *
           std::sqrt(2 * M_PI * C * C /
                     (M * bs.root() * cav.omega(q) * cav.L())) *
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
   * @param q photon momentum
   * @param theta_q angle of photon momentum w.r.t \f$v_s\f$
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
   * where \f$V\f$ is a real matrix describing the polarization vectors (c.f.
   * \ref Cavity::polarizations()) so the final expresion is
   *
   * \f[
   * \tilde{\Pi} = \frac{4\pi c^2}{\omega_q}\frac{2}{L} V^T \Pi(q) V
   * \f]
   *
   * @note that previously we absorbed a factor of \f$e\f$ into the photon
   * fields when defining \f$\pi\f$ so we must restore a factor of \f$e^2 =
   * c\alpha\f$. Thankfully in our units, \f$e=1\f$.
   */
  Matrix2d photon_se_mode(double omega, double q, double theta_q) const
  {
    Matrix2d V = cav.polarizations(q, theta_q);

    auto Piq =
      (Matrix2d() << _photon_se_or_deriv(omega, q, theta_q, 0, 0, false),
       _photon_se_or_deriv(omega, q, theta_q, 0, 1, false),
       _photon_se_or_deriv(omega, q, theta_q, 1, 0, false),
       _photon_se_or_deriv(omega, q, theta_q, 1, 1, false))
        .finished();

    auto Pi_qT = (Matrix2d() << _photon_se_or_deriv(
                    -omega, q, theta_q + M_PI, 0, 0, false),
                  _photon_se_or_deriv(-omega, q, theta_q + M_PI, 0, 1, false),
                  _photon_se_or_deriv(-omega, q, theta_q + M_PI, 1, 0, false),
                  _photon_se_or_deriv(-omega, q, theta_q + M_PI, 1, 1, false))
                   .finished()
                   .transpose();

    return 4 * M_PI * C * C / (cav.L() * cav.omega(q)) * V.transpose() *
           (Piq + Pi_qT) * V;
  }

  /**
   * @brief The derivative of the mode self energy w.r.t frequency
   *
   * @param omega photon frequency
   * @param q photon momentum
   * @param theta_q angle of photon momentum w.r.t \f$v_s\f$
   * @return Matrix2d
   */
  Matrix2d d_photon_se_mode(double omega, double q, double theta_q) const
  {
    Matrix2d V = cav.polarizations(q, theta_q);
    auto Piq =
      (Matrix2d() << _photon_se_or_deriv(omega, q, theta_q, 0, 0, true),
       _photon_se_or_deriv(omega, q, theta_q, 0, 1, true),
       _photon_se_or_deriv(omega, q, theta_q, 1, 0, true),
       _photon_se_or_deriv(omega, q, theta_q, 1, 1, true))
        .finished();

    auto Pi_qT =
      (Matrix2d() << _photon_se_or_deriv(-omega, q, theta_q + M_PI, 0, 0, true),
       _photon_se_or_deriv(-omega, q, theta_q + M_PI, 0, 1, true),
       _photon_se_or_deriv(-omega, q, theta_q + M_PI, 1, 0, true),
       _photon_se_or_deriv(-omega, q, theta_q + M_PI, 1, 1, true))
        .finished()
        .transpose();

    return 4 * M_PI * C * C / (cav.L() * cav.omega(q)) * V.transpose() *
           (Piq + Pi_qT) * V;
  }

  Matrix2d Z(double q, double theta_q, double paraX) const
  {
    return Matrix2d::Identity() -
           paraX * paraX * d_photon_se_mode(cav.omega0, q, theta_q);
  }

  /**
   * @brief Calculate the wavefunction renormalization due to the self energy
   *
   * @param q momentum
   * @param theta_q angle
   * @param paraX enhancement of paramagnetic couping
   * @return auto LLT decomposition of Z
   */
  auto wf_renorm(double q, double theta_q, double paraX) const
  {
    return Z(q, theta_q, paraX).llt();
  }

  /**
   * @}
   */
};
