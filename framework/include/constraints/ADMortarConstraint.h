//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MortarConstraintBase.h"

// Forward Declarations
template <ComputeStage>
class ADMortarConstraint;

declareADValidParams(ADMortarConstraint);

template <ComputeStage compute_stage>
class ADMortarConstraint : public MortarConstraintBase
{
public:
  ADMortarConstraint(const InputParameters & parameters);

  void computeResidual(bool has_master) final;

  void computeJacobian(bool has_master) final;

protected:
  /**
   * compute the residual at the quadrature points
   */
  virtual ADReal computeQpResidual(Moose::MortarType mortar_type) = 0;

  /**
   * compute the residual for the specified element type
   */
  void computeResidual(Moose::MortarType mortar_type) override;

  /**
   * compute the residual for the specified element type
   */
  void computeJacobian(Moose::MortarType mortar_type) override;

private:
  /// A dummy object useful for constructing _lambda when not using Lagrange multipliers
  const ADVariableValue _lambda_dummy;

protected:
  /// The LM solution
  const ADVariableValue & _lambda;

  /// The primal solution on the slave side
  const ADVariableValue & _u_slave;

  /// The primal solution on the master side
  const ADVariableValue & _u_master;

  /// The primal solution gradient on the slave side
  const ADVariableGradient & _grad_u_slave;

  /// The primal solution gradient on the master side
  const ADVariableGradient & _grad_u_master;

  usingCoupleableMembers;
};

#define usingMortarConstraintMembers                                                               \
  usingMortarConstraintBaseMembers;                                                                \
  using ADMortarConstraint<compute_stage>::_lambda;                                                \
  using ADMortarConstraint<compute_stage>::_u_slave;                                               \
  using ADMortarConstraint<compute_stage>::_u_master;                                              \
  using ADMortarConstraint<compute_stage>::_grad_u_slave;                                          \
  using ADMortarConstraint<compute_stage>::_grad_u_master
