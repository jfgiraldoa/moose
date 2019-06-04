//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADMortarConstraint.h"

template <ComputeStage>
class EqualValueConstraint;

declareADValidParams(EqualValueConstraint);

/**
 * Constrain the value of a variable to be the same on both sides of an
 * interface.
 */
template <ComputeStage compute_stage>
class EqualValueConstraint : public ADMortarConstraint<compute_stage>
{
public:
  EqualValueConstraint(const InputParameters & parameters);

protected:
  ADReal computeQpResidual(Moose::MortarType mortar_type) final;

  usingMortarConstraintMembers;
};
