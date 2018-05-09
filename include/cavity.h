#pragma once
#include <cmath>

const double C = 1.; // speed of light
class Cavity {
    public:
    double omega0; // Cavity Frequency

    Cavity(double omega0_) : omega0(omega0_) {}

    double omega(double qx, double qy) const {
        return std::hypot(omega0, C*qx, C*qy);
    }

};