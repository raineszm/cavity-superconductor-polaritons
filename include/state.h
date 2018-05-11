#pragma once

#include "meanfield.h"

class State
{
public:
  double T;
  double delta;

  State(double T_, double delta_)
    : T(T_)
    , delta(delta_)
  {}
};

class Solver
{
public:
  virtual State state(double T) const = 0;
  virtual ~Solver() = default;
};

class MeanFieldSolver : public Solver
{

public:
  double Tc;
  System sys;

  MeanFieldSolver(double Tc_, const System& sys_)
    : Tc(Tc_)
    , sys(sys_)
  {}

  State state(double T) const override
  {
    double delta = MeanField(Tc, T, sys).delta;
    return State(T, delta);
  }

  virtual ~MeanFieldSolver() override = default;
};