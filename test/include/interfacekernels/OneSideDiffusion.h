//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceKernel.h"

// Forward Declarations
class OneSideDiffusion;

template <>
InputParameters validParams<OneSideDiffusion>();

/**
 * DG kernel for interfacing diffusion between two variables on adjacent blocks
 */
class OneSideDiffusion : public InterfaceKernel
{
public:
  OneSideDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type);
  virtual Real computeQpJacobian(Moose::DGJacobianType type);

  Real _D;
};

