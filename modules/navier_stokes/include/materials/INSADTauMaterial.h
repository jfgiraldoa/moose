//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "INSADMaterial.h"

template <ComputeStage>
class INSADTauMaterial;

declareADValidParams(INSADTauMaterial);

template <ComputeStage compute_stage>
class INSADTauMaterial : public INSADMaterial<compute_stage>
{
public:
  INSADTauMaterial(const InputParameters & parameters);

protected:
  virtual void computeProperties() override;
  virtual void computeQpProperties() override;
  void computeHMax();

  const Real _alpha;
  ADMaterialProperty(Real) & _tau;

  ADReal _hmax;

  usingINSMaterialMembers;
};

