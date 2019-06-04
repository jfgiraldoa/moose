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
class MortarConstraint;

template <>
InputParameters validParams<MortarConstraint>();

class MortarConstraint : public MortarConstraintBase
{
public:
  MortarConstraint(const InputParameters & parameters);

  // Using declarations necessary to pull in computeResidual with different parameter list and avoid
  // hidden method warning
  using MortarConstraintBase::computeResidual;

  // Using declarations necessary to pull in computeJacobian with different parameter list and avoid
  // hidden method warning
  using MortarConstraintBase::computeJacobian;

protected:
  /**
   * compute the residual at the quadrature points
   */
  virtual Real computeQpResidual(Moose::MortarType mortar_type) = 0;

  /**
   * compute the jacobian at the quadrature points
   */
  virtual Real computeQpJacobian(Moose::ConstraintJacobianType jacobian_type,
                                 unsigned int jvar) = 0;

  void computeResidual(Moose::MortarType mortar_type) override;

  void computeJacobian(Moose::MortarType mortar_type) override;

private:
  /// A dummy object useful for constructing _lambda when not using Lagrange multipliers
  const VariableValue _lambda_dummy;

protected:
  /// The LM solution
  const VariableValue & _lambda;

  /// The primal solution on the slave side
  const VariableValue & _u_slave;

  /// The primal solution on the master side
  const VariableValue & _u_master;

  /// The primal solution gradient on the slave side
  const VariableGradient & _grad_u_slave;

  /// The primal solution gradient on the master side
  const VariableGradient & _grad_u_master;

  /// The current shape functions
  const VariablePhiValue * _phi;

  /// The current shape function gradients
  const VariablePhiGradient * _grad_phi;

  /// The current shape functions for vector variables
  const VectorVariablePhiValue * _vector_phi;

  /// The current shape function gradients for vector variables
  const VectorVariablePhiGradient * _vector_grad_phi;
};
