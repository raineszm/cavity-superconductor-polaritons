#include "bs.h"
#include "cavity.h"
#include "coupling.h"
#include "exceptions.h"
#include "polariton.h"
#include "state.h"
#include "system.h"
#include <gsl/gsl_errno.h>
#include <iostream>

int
main(int argc, char** argv)
{
  gsl_set_error_handler(&gsl::error_handler);
  const double TC = 0.003016066400623157;
  const double mu = 0.1;
  const double m = 357699.25713625405;

  const double TREL = 0.4;
  const double VREL = 0.9;

  auto sys0 = System(m, mu, TC, 0, 0);
  auto state0 = State::solve(sys0, TREL * TC);
  auto vc = state0.delta / sys0.kf();

  auto sys = System(m, mu, TC, VREL * vc, M_PI_4);
  auto state = State::solve(sys, TREL * TC);
  auto bs = BS(0.1, state);
  // auto p = Polariton(Coupling(bs, Cavity(0.98 * bs.root())), 10);
  auto p = ModePolariton(Coupling(bs, Cavity(0.98 * bs.root())), 10);

  // auto modes =
  //   p.find_modes(0, 0, 0.8 * bs.root(), 1.2 * bs.root(), 1e-2, 1e-17);
  auto modes = p.bands(0, 0);
  std::cout << modes[0] << ", " << modes[1] << modes[2] << std::endl;

  return 0;
}