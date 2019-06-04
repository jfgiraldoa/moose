//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "VectorInitialCondition.h"
#include "Function.h"

// System includes
#include <string>

// Forward Declarations
class VectorFunctionIC;
class InputParameters;

namespace libMesh
{
class Point;
}

template <>
InputParameters validParams<VectorFunctionIC>();

/**
 * IC that calls vectorValue method of a Function object.
 */
class VectorFunctionIC : public VectorInitialCondition
{
public:
  VectorFunctionIC(const InputParameters & parameters);

  virtual RealVectorValue value(const Point & p) override;

protected:
  /// Optional vectorValue function
  const Function * _function;

  /// Optional component function value
  const Function & _function_x;
  const Function & _function_y;
  const Function & _function_z;
};
