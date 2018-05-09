#pragma once 
#include <cmath>
#include <tuple>
#include <array>

#include <rzmcmt/integrate.h>
#include <rzmcmt/fermi.h>

#include "system.h"
#include "state.h"

using rzmcmt::nf;

class Coupling {
    public:
    const System sys;
    const State state;

    double ImDA_int(double kx, double ky, double omega) {
        double fd = std::cos(kx) - std::cos(ky);
        double x = sys.xi(kx, ky);

        double l = std::hypot(x, state.delta);
        double drift = -sys.As/sys.m*kx; // TODO: fix angle of As

        return fd*(nf(drift - l, state.T) - nf(drift + l, state.T))/((omega*omega - 4*l*l)*l);
    }

    double ImDA(double omega, double T) {
        auto integrand = [this, omega](double v[1], const double k[2]) -> int {
            v[0] = ImDA_int(M_PI * k[0], M_PI * k[1], omega);
            return 0;
        };
        std::array<double, 1> result, err;
        std::tie(result, err) = rzmcmt::integrate<2, 1>(integrand, 0.);

        return state.delta*result[0]*omega/sys.m;
    }

};