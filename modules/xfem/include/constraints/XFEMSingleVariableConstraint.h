//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "ElemElemConstraint.h"
#include "MooseMesh.h"

// Forward Declarations
class XFEMSingleVariableConstraint;

class XFEM;

template <>
InputParameters validParams<XFEMSingleVariableConstraint>();

class XFEMSingleVariableConstraint : public ElemElemConstraint
{
public:
  XFEMSingleVariableConstraint(const InputParameters & parameters);
  virtual ~XFEMSingleVariableConstraint();

protected:
  virtual void reinitConstraintQuadrature(const ElementPairInfo & element_pair_info) override;

  virtual Real computeQpResidual(Moose::DGResidualType type) override;

  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

  /// Vector normal to the internal interface
  Point _interface_normal;

  /// Stabilization parameter in Nitsche's formulation
  Real _alpha;

  /// Vector normal to the internal interface
  Real _jump;

  /// Vector normal to the internal interface
  Real _jump_flux;

  /// Use penalty formulation
  bool _use_penalty;

  /// Pointer to the XFEM controller object
  std::shared_ptr<XFEM> _xfem;
};

