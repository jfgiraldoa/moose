//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowMaterialVectorBase.h"
#include "RankTwoTensor.h"

class PorousFlowVolumetricStrain;

template <>
InputParameters validParams<PorousFlowVolumetricStrain>();

/**
 * PorousFlowVolumetricStrain computes volumetric strains, and derivatives thereof
 */
class PorousFlowVolumetricStrain : public PorousFlowMaterialVectorBase
{
public:
  PorousFlowVolumetricStrain(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// If true then the strain rate will include terms that ensure mass is conserved when doing integrals over the displaced mesh
  const bool _consistent;

  /// Number of displacements supplied (1 in 1D, 2 in 2D, 3 in 3D)
  const unsigned int _ndisp;

  /// Displacement variable values at the quad point
  std::vector<const VariableValue *> _disp;

  /// MOOSE variable number of the displacements variables provided
  std::vector<unsigned int> _disp_var_num;

  /// Gradient of the displacements
  std::vector<const VariableGradient *> _grad_disp;

  /// Old value of gradient of the displacements
  std::vector<const VariableGradient *> _grad_disp_old;

  /// The volumetric strain rate at the quadpoints
  MaterialProperty<Real> & _vol_strain_rate_qp;

  /**
   * The derivative of the volumetric strain rate with respect to the porous flow variables.
   * Since the volumetric strain rate depends on derivatives of the displacement variables,
   * this should be multiplied by _grad_phi in kernels
   */
  MaterialProperty<std::vector<RealGradient>> & _dvol_strain_rate_qp_dvar;

  /// The total volumetric strain at the quadpoints
  MaterialProperty<Real> & _vol_total_strain_qp;

  /**
   * The derivative of the total volumetric strain with respect to the porous flow variables.
   * Since the total volumetric strain depends on derivatives of the displacement variables,
   * this should be multiplied by _grad_phi in kernels
   */
  MaterialProperty<std::vector<RealGradient>> & _dvol_total_strain_qp_dvar;
};

