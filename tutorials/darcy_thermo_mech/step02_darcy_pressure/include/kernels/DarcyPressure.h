//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// Including the "Diffusion" Kernel here so we can extend it
#include "Diffusion.h"

class DarcyPressure;

template <>
InputParameters validParams<DarcyPressure>();

/**
 * Computes the residual contribution: K / mu * grad_u * grad_phi.
 *
 * We are inheriting from Diffusion instead of from Kernel because
 * the grad_u * grad_phi is already coded in there and all we
 * need to do is specialize that calculation by multiplying by K / mu.
 */
class DarcyPressure : public Diffusion
{
public:
  DarcyPressure(const InputParameters & parameters);

protected:
  /**
   * Kernels _must_ override computeQpResidual()
   */
  virtual Real computeQpResidual() override;

  /**
   * This is optional (but recommended!)
   */
  virtual Real computeQpJacobian() override;

  /// Will be set from the input file
  Real _permeability;
  Real _viscosity;
};

