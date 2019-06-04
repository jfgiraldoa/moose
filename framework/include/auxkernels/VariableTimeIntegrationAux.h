//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

// forward declarations
class VariableTimeIntegrationAux;

template <>
InputParameters validParams<VariableTimeIntegrationAux>();

/**
 * An AuxKernel that can be used to integrate a field variable in time
 * using a variety of different integration methods.  The result is
 * stored in another field variable.
 */
class VariableTimeIntegrationAux : public AuxKernel
{
public:
  VariableTimeIntegrationAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  Real getIntegralValue();

  std::vector<const VariableValue *> _coupled_vars;
  Real _coef;
  unsigned int _order;
  std::vector<Real> _integration_coef;
};

