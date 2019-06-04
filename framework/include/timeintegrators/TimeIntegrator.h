//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "MooseObject.h"
#include "Restartable.h"

// Forward declarations
class TimeIntegrator;
class FEProblemBase;
class SystemBase;
class NonlinearSystemBase;

namespace libMesh
{
template <typename T>
class NumericVector;
class NonlinearImplicitSystem;
} // namespace libMesh

template <>
InputParameters validParams<TimeIntegrator>();

/**
 * Base class for time integrators
 *
 * Time integrators fulfill two functions:
 * 1) computing u_dot vector (used for computing time derivatives in kernels) and its derivative
 * 2) combining the residual vectors into the final one
 *
 * Capability (1) is used by both NonlinearSystem and AuxiliarySystem, while (2) can be obviously
 * used
 * only by NonlinearSystem (AuxiliarySystem does not produce residual).
 */
class TimeIntegrator : public MooseObject, public Restartable
{
public:
  TimeIntegrator(const InputParameters & parameters);

  /**
   * Called to setup datastructures.
   *
   * WILL be called during recover/restart
   *
   * Should NOT do any computation for initial values
   * use init() for that
   *
   * Note: this doesn't inherit this from SetupInterface because
   * I'm not sure that we need all the other setup functions for
   * TimeIntegrator
   */
  virtual void initialSetup() {}

  /**
   * Called _only_ before the very first timestep (t_step = 0)
   * Never called again (not even during recover/restart)
   */
  virtual void init() {}
  virtual void preSolve() {}
  virtual void preStep() {}

  /**
   * Solves the time step and sets the number of nonlinear and linear iterations.
   */
  virtual void solve();

  /**
   * Callback to the TimeIntegrator called immediately after the
   * residuals are computed in NonlinearSystem::computeResidual().
   * The residual vector which is passed in to this function should
   * be filled in by the user with the _Re_time and _Re_non_time
   * vectors in a way that makes sense for the particular
   * TimeIntegration method.
   */
  virtual void postResidual(NumericVector<Number> & /*residual*/) {}

  /**
   * Callback to the TimeIntegrator called immediately after
   * TimeIntegrator::solve() (so the name does make sense!).  See
   * e.g. CrankNicolson for an example of what can be done in the
   * postSolve() callback -- there it is used to move the residual
   * from the "old" timestep forward in time to avoid recomputing it.
   */
  virtual void postSolve() {}

  /**
   * Callback to the TimeIntegrator called at the very end of time step.
   */
  virtual void postStep() {}

  virtual int order() = 0;
  virtual void computeTimeDerivatives() = 0;

  /**
   * method for computing local automatic differentiation time derivatives
   */
  virtual void computeADTimeDerivatives(DualReal & ad_u_dot, const dof_id_type & dof) const = 0;

  /**
   * Gets the total number of nonlinear iterations over all stages of the time step.
   */
  virtual unsigned int getNumNonlinearIterations() const { return _n_nonlinear_iterations; }

  /**
   * Gets the total number of linear iterations over all stages of the time step.
   */
  virtual unsigned int getNumLinearIterations() const { return _n_linear_iterations; }

protected:
  /**
   * Gets the number of nonlinear iterations in the most recent solve.
   */
  unsigned int getNumNonlinearIterationsLastSolve() const;

  /**
   * Gets the number of linear iterations in the most recent solve.
   */
  unsigned int getNumLinearIterationsLastSolve() const;

  FEProblemBase & _fe_problem;
  SystemBase & _sys;
  NonlinearSystemBase & _nl;

  /// Nonlinear implicit system, if applicable; otherwise, nullptr
  const NonlinearImplicitSystem * _nonlinear_implicit_system;

  /// solution vector for \f$ {du^dot}\over{du} \f$
  Real & _du_dot_du;
  /// solution vectors
  const NumericVector<Number> * const & _solution;
  const NumericVector<Number> & _solution_old;
  const NumericVector<Number> & _solution_older;
  //
  int & _t_step;
  //
  Real & _dt;
  Real & _dt_old;

  /// residual vector for time contributions
  NumericVector<Number> & _Re_time;
  /// residual vector for non-time contributions
  NumericVector<Number> & _Re_non_time;

  /// Total number of nonlinear iterations over all stages of the time step
  unsigned int _n_nonlinear_iterations;
  /// Total number of linear iterations over all stages of the time step
  unsigned int _n_linear_iterations;
};

