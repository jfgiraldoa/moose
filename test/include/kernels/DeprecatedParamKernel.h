//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Reaction.h"

// Forward Declarations
class DeprecatedParamKernel;

template <>
InputParameters validParams<DeprecatedParamKernel>();

class DeprecatedParamKernel : public Reaction
{
public:
  DeprecatedParamKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  Real _coef;
};

