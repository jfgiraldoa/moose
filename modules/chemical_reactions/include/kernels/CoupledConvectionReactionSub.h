//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"

class CoupledConvectionReactionSub;

template <>
InputParameters validParams<CoupledConvectionReactionSub>();

/**
 * Convection of primary species in given equilibrium species
 */
class CoupledConvectionReactionSub : public DerivativeMaterialInterface<Kernel>
{
public:
  CoupledConvectionReactionSub(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// Weight of the equilibrium species concentration in the total primary species concentration
  const Real _weight;
  /// Equilibrium constant for the equilibrium species in association form
  const VariableValue & _log_k;
  /// Stoichiometric coefficient of the primary species
  const Real _sto_u;
  /// Stoichiometric coefficients of the coupled primary species
  const std::vector<Real> _sto_v;
  /// Hydraulic conductivity
  const MaterialProperty<Real> & _cond;
  /// Gravity
  const RealVectorValue _gravity;
  /// Fluid density
  const MaterialProperty<Real> & _density;
  /// Pressure gradient
  const VariableGradient & _grad_p;
  /// Pressure variable number
  const unsigned int _pvar;
  /// Coupled primary species variable numbers
  std::vector<unsigned int> _vars;
  /// Coupled primary species concentrations
  std::vector<const VariableValue *> _vals;
  /// Coupled gradients of primary species concentrations
  std::vector<const VariableGradient *> _grad_vals;
  /// Activity coefficient of primary species in the equilibrium species
  const VariableValue & _gamma_u;
  /// Activity coefficients of coupled primary species in the equilibrium species
  std::vector<const VariableValue *> _gamma_v;
  /// Activity coefficient of equilibrium species
  const VariableValue & _gamma_eq;
};

