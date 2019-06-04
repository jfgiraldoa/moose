//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "NodalConstraint.h"

class EqualValueBoundaryConstraint;

template <>
InputParameters validParams<EqualValueBoundaryConstraint>();

class EqualValueBoundaryConstraint : public NodalConstraint
{
public:
  EqualValueBoundaryConstraint(const InputParameters & parameters);

  /**
   * Called on this object when the mesh changes
   */
  virtual void meshChanged() override;

protected:
  /**
   * Update the sets of nodes with constrained DOFs
   */
  void updateConstrainedNodes();

  /**
   * Computes the residual for the current slave node
   */
  virtual Real computeQpResidual(Moose::ConstraintType type) override;

  /**
   * Computes the jacobian for the constraint
   */
  virtual Real computeQpJacobian(Moose::ConstraintJacobianType type) override;

  // Holds the master node id
  unsigned int _master_node_id;
  // Holds the list of slave node ids
  std::vector<unsigned int> _slave_node_ids;
  // Holds the slave node set or side set
  BoundaryName _slave_node_set_id;
  // Penalty if constraint is not satisfied
  Real _penalty;
};

