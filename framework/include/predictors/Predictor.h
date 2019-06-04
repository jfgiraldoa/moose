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
class Predictor;
class FEProblemBase;
class NonlinearSystemBase;

namespace libMesh
{
template <typename T>
class NumericVector;
}

template <>
InputParameters validParams<Predictor>();

/**
 * Base class for predictors.
 */
class Predictor : public MooseObject, public Restartable
{
public:
  Predictor(const InputParameters & parameters);
  virtual ~Predictor();

  virtual int order() { return 0; }
  virtual void timestepSetup();
  virtual bool shouldApply();
  virtual void apply(NumericVector<Number> & sln) = 0;

  virtual NumericVector<Number> & solutionPredictor() { return _solution_predictor; }

protected:
  FEProblemBase & _fe_problem;
  NonlinearSystemBase & _nl;

  int & _t_step;
  Real & _dt;
  Real & _dt_old;
  const NumericVector<Number> & _solution;
  NumericVector<Number> & _solution_old;
  NumericVector<Number> & _solution_older;
  NumericVector<Number> & _solution_predictor;

  /// Amount by which to scale the predicted value.  Must be in [0,1].
  Real _scale;

  /// Times for which the predictor should not be applied
  std::vector<Real> _skip_times;

  /// Old times for which the predictor should not be applied
  std::vector<Real> _skip_times_old;
};

