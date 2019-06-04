//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADTimeDerivative.h"

template <ComputeStage compute_stage>
class ADHeatConductionTimeDerivative;

declareADValidParams(ADHeatConductionTimeDerivative);

template <ComputeStage compute_stage>
class ADHeatConductionTimeDerivative : public ADTimeDerivative<compute_stage>
{
public:
  ADHeatConductionTimeDerivative(const InputParameters & parameters);

protected:
  virtual ADResidual precomputeQpResidual() override;

  /// Specific heat material property
  const ADMaterialProperty(Real) & _specific_heat;

  /// Density material property
  const ADMaterialProperty(Real) & _density;

  usingTimeDerivativeMembers;
};

