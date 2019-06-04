//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementIntegralPostprocessor.h"

class FunctionElementIntegral;

template <>
InputParameters validParams<FunctionElementIntegral>();

/**
 * Integrates a function over elements
 */
class FunctionElementIntegral : public ElementIntegralPostprocessor
{
public:
  FunctionElementIntegral(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  /// Function to integrate
  const Function & _function;
};
