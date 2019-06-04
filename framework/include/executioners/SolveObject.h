//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseObject.h"
#include "PerfGraphInterface.h"

class SolveObject;
class Executioner;
class FEProblemBase;
class DisplacedProblem;
class MooseMesh;
class NonlinearSystemBase;
class AuxiliarySystem;

class SolveObject : public MooseObject, public PerfGraphInterface
{
public:
  SolveObject(Executioner * ex);

  /**
   * Solve routine provided by this object.
   * @return True if solver is converged.
   */
  virtual bool solve() = 0;

  /// Set the inner solve object wrapped by this object.
  virtual void setInnerSolve(SolveObject & solve) { _inner_solve = &solve; }

protected:
  /// Executioner used to construct this
  Executioner & _executioner;
  /// Reference to FEProblem
  FEProblemBase & _problem;
  /// Displaced problem
  std::shared_ptr<DisplacedProblem> _displaced_problem;
  /// Mesh
  MooseMesh & _mesh;
  /// Displaced mesh
  MooseMesh * _displaced_mesh;
  /// Reference to nonlinear system base for faster access
  NonlinearSystemBase & _nl;
  /// Reference to auxiliary system for faster access
  AuxiliarySystem & _aux;
  /// SolveObject wrapped by this solve object
  SolveObject * _inner_solve;
};
