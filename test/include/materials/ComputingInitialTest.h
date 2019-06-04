//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

// Forward Declarations
class ComputingInitialTest;

template <>
InputParameters validParams<ComputingInitialTest>();

/**
 * Empty material for use in simple applications that don't need material properties.
 */
class ComputingInitialTest : public Material
{
public:
  ComputingInitialTest(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

  MaterialProperty<Real> & _thermal_conductivity;
  const MaterialProperty<Real> & _thermal_conductivity_old;
};

