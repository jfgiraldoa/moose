//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "VectorKernel.h"
#include "MaterialProperty.h"

// Forward Declaration
class LinearVectorPoisson;

template <>
InputParameters validParams<LinearVectorPoisson>();

class LinearVectorPoisson : public VectorKernel
{
public:
  LinearVectorPoisson(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  const Function & _x_sln;
  const Function & _y_sln;
  const Real _eps;
};
