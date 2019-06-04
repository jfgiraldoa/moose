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

class VectorDiffusion;

template <>
InputParameters validParams<VectorDiffusion>();

/**
 * This kernel implements the Laplacian operator:
 * $\nabla \vec{u} \cdot \nabla \vec{\phi_i}$
 */
class VectorDiffusion : public VectorKernel
{
public:
  VectorDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
};

