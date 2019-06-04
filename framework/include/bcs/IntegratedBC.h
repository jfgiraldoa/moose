//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "IntegratedBCBase.h"
#include "MooseVariableInterface.h"

// Forward declarations
class IntegratedBC;

template <>
InputParameters validParams<IntegratedBC>();

/**
 * Base class for deriving any boundary condition of a integrated type
 */
class IntegratedBC : public IntegratedBCBase, public MooseVariableInterface<Real>
{
public:
  IntegratedBC(const InputParameters & parameters);

  virtual MooseVariable & variable() override { return _var; }

  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  /**
   * Computes d-ivar-residual / d-jvar...
   */
  virtual void computeJacobianBlock(MooseVariableFEBase & jvar) override;
  /**
   * Deprecated method
   */
  virtual void computeJacobianBlock(unsigned jvar);
  /**
   * Computes jacobian block with respect to a scalar variable
   * @param jvar The number of the scalar variable
   */
  void computeJacobianBlockScalar(unsigned int jvar) override;

protected:
  /**
   * method for computing the residual at quadrature points
   */
  virtual Real computeQpResidual() = 0;

  MooseVariable & _var;

  /// normals at quadrature points
  const MooseArray<Point> & _normals;

  // shape functions

  /// shape function values (in QPs)
  const VariablePhiValue & _phi;
  /// gradients of shape functions (in QPs)
  const VariablePhiGradient & _grad_phi;

  // test functions

  /// test function values (in QPs)
  const VariableTestValue & _test;
  /// gradients of test functions  (in QPs)
  const VariableTestGradient & _grad_test;

  // unknown

  /// the values of the unknown variable this BC is acting on
  const VariableValue & _u;
  /// the gradient of the unknown variable this BC is acting on
  const VariableGradient & _grad_u;
};

