//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"

class Function;
class Gravity;

template <>
InputParameters validParams<Gravity>();

/**
 * Gravity computes the body force (force/volume) given the acceleration of gravity (value) and the
 * density
 */
class Gravity : public Kernel
{
public:
  Gravity(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  const MaterialProperty<Real> & _density;
  const Real _value;
  const Function & _function;

  // _alpha parameter for HHT time integration scheme
  const Real _alpha;
};
