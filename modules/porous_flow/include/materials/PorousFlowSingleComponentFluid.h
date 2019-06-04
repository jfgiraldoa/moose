//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowFluidPropertiesBase.h"
#include "SinglePhaseFluidProperties.h"

class PorousFlowSingleComponentFluid;

template <>
InputParameters validParams<PorousFlowSingleComponentFluid>();

/**
 * General single component fluid material. Provides quadpoint density, viscosity,
 * internal energy, enthalpy and derivatives wrt pressure and temperature
 * for a fluid defined in the FluidProperties module
 */
class PorousFlowSingleComponentFluid : public PorousFlowFluidPropertiesBase
{
public:
  PorousFlowSingleComponentFluid(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// Unit used for porepressure
  const enum class PressureUnitEnum { Pa, MPa } _p_unit;

  /// convert porepressure to Pascals by multiplying by this quantity
  const Real _pressure_to_Pascals;

  /// Unit used for time
  const enum class TimeUnitEnum { seconds, hours, days, years } _time_unit;

  /// convert time to seconds by multiplying by this quantity
  const Real _time_to_seconds;

  /// If true, this Material will compute density and viscosity, and their derivatives
  const bool _compute_rho_mu;

  /// If true, this Material will compute internal energy and its derivatives
  const bool _compute_internal_energy;

  /// If true, this Material will compute enthalpy and its derivatives
  const bool _compute_enthalpy;

  /// Fluid phase density at the qps or nodes
  MaterialProperty<Real> * const _density;

  /// Derivative of fluid density wrt phase pore pressure at the qps or nodes
  MaterialProperty<Real> * const _ddensity_dp;

  /// Derivative of fluid density wrt temperature at the qps or nodes
  MaterialProperty<Real> * const _ddensity_dT;

  /// Fluid phase viscosity at the nodes
  MaterialProperty<Real> * const _viscosity;

  /// Derivative of fluid phase viscosity wrt pressure at the nodes
  MaterialProperty<Real> * const _dviscosity_dp;

  /// Derivative of fluid phase viscosity wrt temperature at the nodes
  MaterialProperty<Real> * const _dviscosity_dT;

  /// Fluid phase internal_energy at the qps or nodes
  MaterialProperty<Real> * const _internal_energy;

  /// Derivative of fluid internal_energy wrt phase pore pressure at the qps or nodes
  MaterialProperty<Real> * const _dinternal_energy_dp;

  /// Derivative of fluid internal_energy wrt temperature at the qps or nodes
  MaterialProperty<Real> * const _dinternal_energy_dT;

  /// Fluid phase enthalpy at the qps or nodes
  MaterialProperty<Real> * const _enthalpy;

  /// Derivative of fluid enthalpy wrt phase pore pressure at the qps or nodes
  MaterialProperty<Real> * const _denthalpy_dp;

  /// Derivative of fluid enthalpy wrt temperature at the qps or nodes
  MaterialProperty<Real> * const _denthalpy_dT;

  /// Fluid properties UserObject
  const SinglePhaseFluidProperties & _fp;
};
