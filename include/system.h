#pragma once

#include "integrate.h"
#include "utils.h"
#include <cmath>

//! The material properties of the system

//! This class is the first instantiated in any code.
class System
{
public:
  //! The electron effective mass in units of m_e
  double m;
  //! The chemical potential
  double mu;
  //! Critical temperature
  double Tc;
  //! The superfluid velocity
  double vs;
  //! The angle that the superfluid velocity makes with the x axis
  double theta_v;

  System(double m_, double mu_, double Tc_, double vs_, double theta_v_)
    : m(m_)
    , mu(mu_)
    , Tc(Tc_)
    , vs(vs_)
    , theta_v(theta_v_)
  {}

  double xi(double kx, double ky) const
  {
    return (kx * kx + ky * ky) / (2 * m) - mu + 0.5 * m * vs * vs;
  }

  //! The "bare" energy plus the superfluid term
  //! \f[\xi_{\mathbf k} = \frac{k^2}{2m} - \mu + \frac{1}{2} m v_s^2\f]
  double xi_k(double k) const
  {
    return k * k / (2 * m) - mu + 0.5 * m * vs * vs;
  }

  //! The Fermi momentum
  double kf() const { return std::sqrt(2 * m * mu); }
  //! The Fermi velocity
  double vf() const { return kf() / m; }
  //! The density of states per spin in 2D \f[\nu=\frac{m}{2\pi}\f]
  double dos() const { return m / (2 * M_PI); }
  /** The density of electrons.
   *
   * For low temperatures we can neglect the distinction between \f$E_f\f$ and
   * \f$\mu\f$. In that case the density is \f[ n = 2 \int
   * \frac{d\mathbf{k}}{2\pi} \Theta(-\xi_\mathbf{k}) \approx 2 \nu
   * \int_{-E_f}^0 d\xi = 2\nu E_f \f]
   */
  constexpr double n() const { return 2 * dos() * mu; }

  double drift(double kx, double ky) const
  {
    return drift_theta(std::hypot(kx, ky), std::atan2(ky, kx));
  }

  //! The Doppler shift term \f$\mathbf{v}_s \cdot \mathbf{k}\f$
  double drift_theta(double k, double theta) const
  {
    return vs * k * std::cos(theta - theta_v);
  }

  /** The gap equation for the s-wave state written in the form
  \f[2\nu \int_0^\infty d\xi \left(\left[
  \int \frac{d\theta}{2\pi} \frac{\tanh\frac{v_s \cdot k + E}{2T} -
  \tanh\frac{v_s \cdot k - E}{2T}}{4E} \right] -
  \frac{\tanh\frac{\xi}{2T_c}}{2\xi}\right)\f]
  */
  double gap_eq(double T, double delta) const
  {
    return 2 * dos() *
           gsl_xi_integrate(
             [this, delta, T](double x, double theta) {
               return gap_eq_int(x, theta, T, delta);
             },
             0.);
  }

  //! The integrand of the gap equation
  double gap_eq_int(double x, double theta, double T, double delta) const
  {

    double l = std::hypot(x, delta);
    double d = drift_theta(kf(), theta);
    return c(l, d, T) - tanh_over(x, Tc) / 2;
  }
};