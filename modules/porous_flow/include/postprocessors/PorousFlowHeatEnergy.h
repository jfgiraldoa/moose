//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementIntegralVariablePostprocessor.h"
#include "PorousFlowDictator.h"

class PorousFlowHeatEnergy;

template <>
InputParameters validParams<PorousFlowHeatEnergy>();

/**
 * Postprocessor produces the sum of heat energy of the porous skeleton and/or fluid components in a
 * region
 */
class PorousFlowHeatEnergy : public ElementIntegralPostprocessor
{
public:
  PorousFlowHeatEnergy(const InputParameters & parameters);

protected:
  virtual Real computeIntegral() override;
  virtual Real computeQpIntegral() override;

  /// PorousFlowDictator UserObject
  const PorousFlowDictator & _dictator;

  /// Number of fluid phases
  const unsigned int _num_phases;

  /// Whether fluid is present
  const bool _fluid_present;

  /// Whether to include the heat energy of the porous skeleton in the calculations
  const bool _include_porous_skeleton;

  /// The phase indices that this Postprocessor is restricted to
  std::vector<unsigned int> _phase_index;

  /// Porosity
  const MaterialProperty<Real> & _porosity;

  /// Nodal rock energy density
  const MaterialProperty<Real> & _rock_energy_nodal;

  /// Nodal fluid density
  const MaterialProperty<std::vector<Real>> * const _fluid_density;

  /// Nodal fluid saturation
  const MaterialProperty<std::vector<Real>> * const _fluid_saturation_nodal;

  /// Internal energy of the phases, evaluated at the nodes
  const MaterialProperty<std::vector<Real>> * const _energy_nodal;

  /// The variable for the corresponding PorousFlowEnergyTimeDerivative Kernel: this provides test functions
  MooseVariable * const _var;
};

