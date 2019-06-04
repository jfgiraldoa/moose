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
#include "LinearInterpolation.h"
#include "DerivativeMaterialInterface.h"

// Forward Declarations
class PiecewiseLinearInterpolationMaterial;

template <>
InputParameters validParams<PiecewiseLinearInterpolationMaterial>();

/**
 * This material uses a LinearInterpolation object to define the dependence
 * of the material's value on a variable.
 */
class PiecewiseLinearInterpolationMaterial : public DerivativeMaterialInterface<Material>
{
public:
  PiecewiseLinearInterpolationMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Name of the property to be computed
  std::string _prop_name;

  /// Value of the coupled variable to be used as the abscissa in the piecewise linear interpolation
  const VariableValue & _coupled_var;

  /// Factor to scale the ordinate values by (default = 1)
  const Real _scale_factor;

  /// Material property to be calculated
  MaterialProperty<Real> & _property;
  /// First derivative of the material property wrt the coupled variable
  MaterialProperty<Real> & _dproperty;

  /// LinearInterpolation object
  std::unique_ptr<LinearInterpolation> _linear_interp;
};

