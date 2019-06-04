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
#include "MooseObject.h"
#include "SetupInterface.h"
#include "FunctionInterface.h"
#include "UserObjectInterface.h"
#include "TransientInterface.h"
#include "GeometricSearchInterface.h"
#include "Restartable.h"
#include "MeshChangedInterface.h"
#include "TaggingInterface.h"

// Forward Declarations
class Assembly;
class Constraint;
template <typename>
class MooseVariableFE;
typedef MooseVariableFE<Real> MooseVariable;
typedef MooseVariableFE<VectorValue<Real>> VectorMooseVariable;
class SubProblem;
class MooseMesh;

template <>
InputParameters validParams<Constraint>();

/**
 * Base class for all Constraint types
 */
class Constraint : public MooseObject,
                   public SetupInterface,
                   public FunctionInterface,
                   public UserObjectInterface,
                   public TransientInterface,
                   protected GeometricSearchInterface,
                   public Restartable,
                   public MeshChangedInterface,
                   public TaggingInterface
{
public:
  Constraint(const InputParameters & parameters);
  virtual ~Constraint();

  /**
   * Subproblem this constraint is part of
   * @return The reference to the subproblem
   */
  SubProblem & subProblem() { return _subproblem; }

  virtual bool addCouplingEntriesToJacobian() { return true; }
  virtual void subdomainSetup() override final
  {
    mooseError("subdomain setup for constraints is not implemented");
  }

  virtual void residualEnd() {}

protected:
  SystemBase & _sys;

  THREAD_ID _tid;

  Assembly & _assembly;
  MooseMesh & _mesh;

  unsigned int _i, _j;
  unsigned int _qp;
};

#define usingConstraintMembers                                                                     \
  usingMooseObjectMembers;                                                                         \
  usingTaggingInterfaceMembers;                                                                    \
  using Constraint::_i;                                                                            \
  using Constraint::_qp;                                                                           \
  using Constraint::_tid
