//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementAverageValue.h"

// Forward Declarations
class ElementAverageSecondTimeDerivative;

template <>
InputParameters validParams<ElementAverageSecondTimeDerivative>();

/**
 * This postprocessor computes a volume integral of the second time derivative of a given variable.
 */
class ElementAverageSecondTimeDerivative : public ElementAverageValue
{
public:
  ElementAverageSecondTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  /// Holds the solution second derivative wrt time at the current quadrature points
  const VariableValue & _u_dotdot;
};

