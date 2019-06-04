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

// Forward Declarations
class NormalizationAux;

template <>
InputParameters validParams<NormalizationAux>();

/**
 * This auxiliary kernel normalizes a variable based on a postprocessor.
 * Typically this postprocessor is a norm of the variable to be normalized.
 * The option to shift the value is also provided.
 */
class NormalizationAux : public AuxKernel
{
public:
  NormalizationAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const VariableValue & _src;
  const Real * _pp_on_source;
  const Real * _shift;
  Real _normal_factor;
};

