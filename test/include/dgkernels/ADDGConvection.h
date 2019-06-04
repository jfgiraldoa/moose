//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADDGKernel.h"

template <ComputeStage>
class ADDGConvection;

declareADValidParams(ADDGConvection);

template <ComputeStage compute_stage>
class ADDGConvection : public ADDGKernel<compute_stage>
{
public:
  ADDGConvection(const InputParameters & parameters);

protected:
  ADReal computeQpResidual(Moose::DGResidualType type) override;

  RealVectorValue _velocity;

  usingDGKernelMembers;
};
