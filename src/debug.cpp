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
  const double TC = 0.8186465944548569;
  const double mu = 5325.511530243174;
  const double m = 1.6;

  const double TREL = 0.4;
  const double VREL = 0.6;

  auto sys0 = System(m, mu, TC, 0, 0);
  auto state0 = State::solve(sys0, TREL * TC);
  auto vc = state0.delta / sys0.kf();

  auto sys = System(m, mu, TC, VREL * vc, 0);
  auto state = State::solve(sys, TREL * TC);
  auto bs = BS(state.delta / 10, state);
  auto p = Polariton(bs, Cavity(1.04 * bs.root()), Coupling(state));

  auto modes = p.find_modes(0, 0);
  std::cout << modes[0] << ", " << modes[1] << modes[2] << std::endl;

  return 0;
}