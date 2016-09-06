/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowDispersiveFlux.h"

template<>
InputParameters validParams<PorousFlowDispersiveFlux>()
{
  InputParameters params = validParams<Kernel>();
  params.addParam<unsigned int>("fluid_component", 0, "The index corresponding to the fluid component for this kernel");
  params.addRequiredParam<UserObjectName>("PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names");
  params.addRequiredParam<std::vector<Real> >("disp_long", "Vector of longitudinal dispersion coefficients for each phase");
  params.addRequiredParam<std::vector<Real> >("disp_trans", "Vector of transverse dispersion coefficients for each phase");
  params.addRequiredParam<RealVectorValue>("gravity", "Gravitational acceleration vector downwards (m/s^2)");
  params.addClassDescription("Dispersive and diffusive flux of the component given by fluid_component in all phases");
  return params;
}

PorousFlowDispersiveFlux::PorousFlowDispersiveFlux(const InputParameters & parameters) :
    Kernel(parameters),

    _fluid_density_qp(getMaterialProperty<std::vector<Real> >("PorousFlow_fluid_phase_density_qp")),
    _dfluid_density_qp_dvar(getMaterialProperty<std::vector<std::vector<Real> > >("dPorousFlow_fluid_phase_density_qp_dvar")),
    _grad_mass_frac(getMaterialProperty<std::vector<std::vector<RealGradient> > >("PorousFlow_grad_mass_frac")),
    _dmass_frac_dvar(getMaterialProperty<std::vector<std::vector<std::vector<Real> > > >("dPorousFlow_mass_frac_dvar")),
    _porosity_qp(getMaterialProperty<Real>("PorousFlow_porosity_qp")),
    _dporosity_qp_dvar(getMaterialProperty<std::vector<Real> >("dPorousFlow_porosity_qp_dvar")),
    _tortuosity(getMaterialProperty<std::vector<Real> >("PorousFlow_tortuosity")),
    _dtortuosity_dvar(getMaterialProperty<std::vector<std::vector<Real> > >("dPorousFlow_tortuosity_dvar")),
    _diffusion_coeff(getMaterialProperty<std::vector<std::vector<Real> > >("PorousFlow_diffusion_coeff")),
    _ddiffusion_coeff_dvar(getMaterialProperty<std::vector<std::vector<std::vector<Real> > > >("dPorousFlow_diffusion_coeff_dvar")),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _fluid_component(getParam<unsigned int>("fluid_component")),
    _num_phases(_dictator.numPhases()),
    _identity_tensor(RankTwoTensor::initIdentity),
    _relative_permeability(getMaterialProperty<std::vector<Real> >("PorousFlow_relative_permeability")),
    _drelative_permeability_dvar(getMaterialProperty<std::vector<std::vector<Real> > >("dPorousFlow_relative_permeability_dvar")),
    _fluid_viscosity(getMaterialProperty<std::vector<Real> >("PorousFlow_viscosity")),
    _dfluid_viscosity_dvar(getMaterialProperty<std::vector<std::vector<Real> > >("dPorousFlow_viscosity_dvar")),
    _permeability(getMaterialProperty<RealTensorValue>("PorousFlow_permeability_qp")),
    _dpermeability_dvar(getMaterialProperty<std::vector<RealTensorValue> >("dPorousFlow_permeability_qp_dvar")),
    _grad_p(getMaterialProperty<std::vector<RealGradient> >("PorousFlow_grad_porepressure_qp")),
    _dgrad_p_dgrad_var(getMaterialProperty<std::vector<std::vector<Real> > >("dPorousFlow_grad_porepressure_qp_dgradvar")),
    _dgrad_p_dvar(getMaterialProperty<std::vector<std::vector<RealGradient> > >("dPorousFlow_grad_porepressure_qp_dvar")),
    _gravity(getParam<RealVectorValue>("gravity")),
    _disp_long(getParam<std::vector<Real> >("disp_long")),
    _disp_trans(getParam<std::vector<Real> >("disp_trans"))
  {
  // Check that sufficient values of the dispersion coefficients have been entered
  if (_disp_long.size() != _num_phases)
    mooseError("The number of longitudinal dispersion coefficients disp_long in " << _name << " is not equal to the number of phases");

  if (_disp_trans.size() != _num_phases)
    mooseError("The number of transverse dispersion coefficients disp_trans in " << _name << " is not equal to the number of phases");
}

Real
PorousFlowDispersiveFlux::computeQpResidual()
{
  RealVectorValue flux = 0.0;
  RealVectorValue velocity;
  Real velocity_abs;
  RankTwoTensor v2;
  RankTwoTensor dispersion;
  dispersion.zero();
  Real diffusion;

  for (unsigned int ph = 0; ph < _num_phases; ++ph)
  {
    // Diffusive component
    diffusion = _porosity_qp[_qp] * _tortuosity[_qp][ph] * _diffusion_coeff[_qp][ph][_fluid_component];

    // Calculate Darcy velocity
    velocity = (_permeability[_qp] * (_grad_p[_qp][ph] - _fluid_density_qp[_qp][ph] *
      _gravity) * _relative_permeability[_qp][ph] / _fluid_viscosity[_qp][ph]);
    velocity_abs = std::sqrt(velocity * velocity);

    if (velocity_abs > 0.0)
    {
      v2.vectorOuterProduct(velocity, velocity);

      // Add longitudinal dispersion to diffusive component
      diffusion += _disp_trans[ph] * velocity_abs;
      dispersion = (_disp_long[ph] - _disp_trans[ph]) * v2 / velocity_abs;
    }

    flux += _fluid_density_qp[_qp][ph] * (diffusion * _identity_tensor + dispersion) * _grad_mass_frac[_qp][ph][_fluid_component];
  }
  return _grad_test[_i][_qp] * flux;
}

Real
PorousFlowDispersiveFlux::computeQpJacobian()
{
  return computeQpJac(_var.number());
}

Real
PorousFlowDispersiveFlux::computeQpOffDiagJacobian(unsigned int jvar)
{
  return computeQpJac(jvar);
}

Real
PorousFlowDispersiveFlux::computeQpJac(unsigned int jvar) const
{
  // If the variable is not a valid PorousFlow variable, set the Jacobian to 0
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;

  const unsigned int pvar = _dictator.porousFlowVariableNum(jvar);

  RealVectorValue velocity;
  Real velocity_abs;
  RankTwoTensor v2;
  RankTwoTensor dispersion;
  dispersion.zero();
  Real diffusion;
  RealVectorValue flux = 0.0;
  RealVectorValue dflux = 0.0;

  for (unsigned int ph = 0; ph < _num_phases; ++ph)
  {
    // Diffusive component
    diffusion = _porosity_qp[_qp] * _tortuosity[_qp][ph] * _diffusion_coeff[_qp][ph][_fluid_component];

    // Calculate Darcy velocity
    velocity = (_permeability[_qp] * (_grad_p[_qp][ph] - _fluid_density_qp[_qp][ph] *
      _gravity) * _relative_permeability[_qp][ph] / _fluid_viscosity[_qp][ph]);
    velocity_abs = std::sqrt(velocity * velocity);

    if (velocity_abs > 0.0)
    {
      v2.vectorOuterProduct(velocity, velocity);

      // Add longitudinal dispersion to diffusive component
      diffusion += _disp_trans[ph] * velocity_abs;
      dispersion = (_disp_long[ph] - _disp_trans[ph]) * v2 / velocity_abs;
    }

    // Derivative of Darcy velocity
    RealVectorValue dvelocity = _dpermeability_dvar[_qp][pvar] * (_grad_p[_qp][ph] - _fluid_density_qp[_qp][ph]*_gravity);
    dvelocity += _permeability[_qp] * (_grad_phi[_j][_qp] * _dgrad_p_dgrad_var[_qp][ph][pvar] - _phi[_j][_qp] * _dfluid_density_qp_dvar[_qp][ph][pvar] * _gravity);
    dvelocity += _permeability[_qp] * (_dgrad_p_dvar[_qp][ph][pvar] * _phi[_j][_qp]);

    Real dvelocity_abs = 0.0;
    if (velocity_abs > 0.0)
      dvelocity_abs = velocity * dvelocity / velocity_abs;

    // Derivative of diffusion term (note: dispersivity is assumed constant)
    Real ddiffusion = _phi[_j][_qp] * _dporosity_qp_dvar[_qp][pvar] * _tortuosity[_qp][ph] * _diffusion_coeff[_qp][ph][_fluid_component];
    ddiffusion += _phi[_j][_qp] * _porosity_qp[_qp] * _dtortuosity_dvar[_qp][ph][pvar] * _diffusion_coeff[_qp][ph][_fluid_component];
    ddiffusion += _phi[_j][_qp] * _porosity_qp[_qp] * _tortuosity[_qp][ph] * _ddiffusion_coeff_dvar[_qp][ph][_fluid_component][pvar];
    ddiffusion += _disp_trans[ph] * dvelocity_abs;

    // Derivative of dispersion term (note: dispersivity is assumed constant)
    RankTwoTensor ddispersion;
    ddispersion.zero();
    if (velocity_abs > 0.0)
    {
      RankTwoTensor dv2a, dv2b;
      dv2a.vectorOuterProduct(velocity, dvelocity);
      dv2b.vectorOuterProduct(dvelocity, velocity);
      ddispersion = (_disp_long[ph] - _disp_trans[ph]) * (dv2a + dv2b) / velocity_abs;
      ddispersion -= (_disp_long[ph] - _disp_trans[ph]) * v2 * dvelocity_abs / velocity_abs / velocity_abs;
    }

    dflux += _phi[_j][_qp] * _dfluid_density_qp_dvar[_qp][ph][pvar] * (diffusion * _identity_tensor + dispersion) * _grad_mass_frac[_qp][ph][_fluid_component];
    dflux += _fluid_density_qp[_qp][ph] * (ddiffusion * _identity_tensor + ddispersion) * _grad_mass_frac[_qp][ph][_fluid_component];
    dflux += _fluid_density_qp[_qp][ph] * (diffusion * _identity_tensor + dispersion) * _dmass_frac_dvar[_qp][ph][_fluid_component][pvar] * _grad_phi[_j][_qp];
  }

  return _grad_test[_i][_qp] * dflux;
}