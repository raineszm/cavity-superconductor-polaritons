//! \file cavity.h
//!  We use Gaussian Hartree Units
#pragma once
#include <cmath>

//! The fine structure constant \f$\alpha\f$
const double ALPHA = 7.2973525664e-3;

//! The speed of light
const double C = 1 / ALPHA;

//! The paramagnetic coupling strength
const double GPAR = std::sqrt(ALPHA / C);

//! The physics of the cavity modes
class Cavity
{
public:
  //! Cavity frequency
  double omega0;

  explicit Cavity(double omega0_)
    : omega0(omega0_)
  {}

  //! The dispersion of the cavity
  //! \f[\omega^2 = \omega_0^2 + c^2 q^2\f]
  double omega(double qx, double qy) const
  {
    return std::sqrt(omega0 * omega0 + C * C * (qx * qx + qy * qy));
  }

  // TODO: handle polarization stuff
  //! The Lagrangian of the photon modes

  //! \f[\mathcal{L} = \omega^2 - \omega_0^2 - c^2q^2\f]
  //! gives rise to action \f$S = \frac{1}{8\pi c^2}\int A_q\mathcal{L}
  //! A_{-q}\f$
  double action(double omega, double qx, double qy) const
  {
    return (omega * omega - this->omega(qx, qy)) / (8 * M_PI * C * C);
  }
};