//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "TensorMechanicsHardeningModel.h"

class TensorMechanicsHardeningConstant;

template <>
InputParameters validParams<TensorMechanicsHardeningConstant>();

/**
 * No hardening - the parameter assumes the value _val
 * for all internal parameters
 */
class TensorMechanicsHardeningConstant : public TensorMechanicsHardeningModel
{
public:
  TensorMechanicsHardeningConstant(const InputParameters & parameters);

  virtual Real value(Real intnl) const override;

  virtual Real derivative(Real intnl) const override;

  virtual std::string modelName() const override;

private:
  /// The value that the parameter will take
  Real _val;
};

