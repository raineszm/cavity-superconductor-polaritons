/**
 * @brief The electronic system
 *
 * @file system.h
 * @author Zach Raines
 * @date 2018-06-12
 */
#pragma once

#include "integrate.h"
#include "utils.h"
#include <array>
#include <cmath>
#include <functional>

/**
 * @brief The material properties of the system
 *
 * This class is the first instantiated in any code.
 */
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
  double theta_s;

  /**
   * @brief Construct a new System object
   *
   * @param m_ #m
   * @param mu_ #mu
   * @param Tc_ #Tc
   * @param vs_ #vs
   * @param theta_s_ #theta_s
   */
  System(double m_, double mu_, double Tc_, double vs_, double theta_s_)
    : m(m_)
    , mu(mu_)
    , Tc(Tc_)
    , vs(vs_)
    , theta_s(theta_s_)
  {}

  bool operator==(const System& rhs) const
  {
    return m == rhs.m && mu == rhs.mu && Tc == rhs.Tc && vs == rhs.vs &&
           theta_s == rhs.theta_s;
  }

  /**
   * @brief Bare qp energy in cartesian coordinates
   *
   * @param kx
   * @param ky
   * @return double
   */
  double xi(double kx, double ky) const { return xi_k(std::hypot(kx, ky)); }

  /** @name Dispersion
   * @{
   */

  /**
   * @brief The "bare" energy plus the superfluid term
   *
   * @param k the magnitude of the electronic momentum
   * @return double
   *
   * \f[\xi_{\mathbf k} = \frac{k^2}{2m} - \mu + \frac{1}{2} m v_s^2\f]
   */
  double xi_k(double k) const
  {
    return k * k / (2 * m) - mu + 0.5 * m * vs * vs;
  }

  //! The Fermi momentum
  double kf() const { return std::sqrt(2 * m * mu); }
  //! The Fermi velocity
  double vf() const { return kf() / m; }
  //! The density of states per spin in 2D \f[\nu=\frac{m}{2\pi}\f]
  constexpr double dos() const { return m / (2 * M_PI); }

  /** @} */

  /**
   * @brief The density of electrons.
   *
   * @return constexpr double
   *
   * For low temperatures we can neglect the distinction between \f$E_f\f$ and
   * \f$\mu\f$. In that case the density is \f[ n = 2 \int
   * \frac{d\mathbf{k}}{2\pi} \Theta(-\xi_\mathbf{k}) \approx 2 \nu
   * \int_{-E_f}^0 d\xi = 2\nu E_f \f]
   */
  constexpr double n() const { return 2 * dos() * mu; }

  /**
   * @name Doppler shift
   *
   * @{
   */

  /**
   * @brief The Doppler shift term \f$\mathbf{v}_s \cdot \mathbf{k}\f$
   *
   * @param k the magnitude of the electron momentum
   * @param theta the angle the momentum makes with the x axis
   * @return double
   */
  double doppler(double k, double theta) const
  {
    return vs * k * std::cos(theta);
  }

  /**
   * @}
   */

  /** @name Gap Equation
   * @{
   */

  /**
   * @brief The gap equation for the s-wave state written in the form
   *
   * @param T temperature
   * @param delta the gap
   * @return double
   *
   * \f[2\nu \int_0^\infty d\xi \left(\left[
   * \int \frac{d\theta}{2\pi} \frac{\tanh\frac{v_s \cdot k + E}{2T} -
   * \tanh\frac{v_s \cdot k - E}{2T}}{4E} \right] -
   * \frac{\tanh\frac{\xi}{2T_c}}{2\xi}\right)\f]
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

  /**
   * @brief The integrand appearing in the gap equation
   *
   * @param x the quasiparticle energy
   * @param theta angle on the Fermi surface
   * @param T temperature
   * @param delta gap size
   * @return double
   */
  double gap_eq_int(double x, double theta, double T, double delta) const
  {

    double l = std::hypot(x, delta);
    double d = doppler(kf(), theta);
    return c(l, d, T) - tanh_over(x, Tc) / 2;
  }

  /** @} */

  /**
   * @brief Select the \f$i\f$ components of \f$\mathbf{v}_s\f$
   *
   * @param i component of vector
   * @return double
   * @see v_comp()
   *
   * Here we operate in the basis where the \f$x\f$ axis is along
   * \f$\mathbf{v}_s\f$.
   */
  double vs_comp(int i) const
  {
    assert(i < 2 and i >= 0);
    if (i == 0) {
      return vs;
    } else {
      return 0.;
    }
  }

  /**
   * @brief Select the \f$i\f$ components of \f$\mathbf{v}\f$
   *
   * @param kx x component of momentum
   * @param ky y component of momentum
   * @param i component of vector
   * @return double
   * @see vs_comp()
   */
  double v_comp(double k, double theta, int i) const
  {
    assert(i < 2 and i >= 0);
    return gsl::polar_to_rect(k / m, theta)[i];
  }

  using pickle_type = std::array<double, 5>;

  pickle_type pickle() const { return { { m, mu, Tc, vs, theta_s } }; }

  static inline System unpickle(const pickle_type& a)
  {
    return std::make_from_tuple<System, decltype(a)>(a);
  }
};

//! \cond
namespace std {
template<>
struct hash<System>
{
  typedef System argument_type;
  typedef std::size_t result_type;
  result_type operator()(argument_type const& s) const noexcept
  {
    result_type const hm = std::hash<double>{}(s.m);
    result_type const hmu = std::hash<double>{}(s.mu);
    result_type const hTc = std::hash<double>{}(s.Tc);
    result_type const hvs = std::hash<double>{}(s.vs);
    result_type const htheta_s = std::hash<double>{}(s.theta_s);
    return hm ^ (hmu << 1) ^ (hTc << 2) ^ (hvs << 3) ^ (htheta_s << 4);
  }
};
}
//! \endcond
