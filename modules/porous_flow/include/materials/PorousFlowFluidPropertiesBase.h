//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowMaterialBase.h"
#include "PorousFlowDictator.h"

class PorousFlowFluidPropertiesBase;

template <>
InputParameters validParams<PorousFlowFluidPropertiesBase>();

/**
 * Base class for fluid properties materials. All PorousFlow fluid
 * materials must override computeQpProperties()
 */
class PorousFlowFluidPropertiesBase : public PorousFlowMaterialBase
{
public:
  PorousFlowFluidPropertiesBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Pore pressure at the nodes or quadpoints
  const MaterialProperty<std::vector<Real>> & _porepressure;

  /// Fluid temperature at the nodes or quadpoints
  const MaterialProperty<Real> & _temperature;

  /// Conversion from degrees Celsius to degrees Kelvin
  const Real _t_c2k;

  /// Universal gas constant
  const Real _R;
};

