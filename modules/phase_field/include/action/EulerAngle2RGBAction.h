//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InputParameters.h"
#include "Action.h"

/**
 * Automatically generates all variables, Kernels, and Materials to ensure the
 * correct derivatives of the elastic free energy in a non-split Cahn-Hilliard
 * simulation are assembled.
 */
class EulerAngle2RGBAction : public Action
{
public:
  EulerAngle2RGBAction(const InputParameters & params);

  virtual void act();

private:
  const std::string _var_name_base;
};

template <>
InputParameters validParams<EulerAngle2RGBAction>();

