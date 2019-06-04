//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "IntegratedBC.h"

// Forward Declarations
class FunctionGradientNeumannBC;
class Function;

template <>
InputParameters validParams<FunctionGradientNeumannBC>();

/**
 * Boundary condition of a Neumann style whose value is computed by a user-defined function
 */
class FunctionGradientNeumannBC : public IntegratedBC
{
public:
  FunctionGradientNeumannBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  /// The function being used for setting the value
  const Function & _exact_solution;
  const Real _coeff;
};
