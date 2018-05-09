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
  virtual double Tc() const = 0;
  virtual ~Solver() = default;
};

class MeanFieldSolver : public Solver
{

public:
  double g;
  System sys;

  MeanFieldSolver(double g_, const System& sys_)
    : g(g_)
    , sys(sys_)
  {}

  State state(double T) const override
  {
    double delta = MeanField(g, T, sys).delta;
    return State(T, delta);
  }

  double Tc() const override
  {
    return MeanField::Tc(g, sys);
  }

  virtual ~MeanFieldSolver() override = default;
};