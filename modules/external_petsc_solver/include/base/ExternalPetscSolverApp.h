//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseApp.h"

class ExternalPetscSolverApp;

template <>
InputParameters validParams<ExternalPetscSolverApp>();

/**
 * This is a demo used to demonstrate how to couple an external app
 * through the MOOSE wrapper APP.
 * We are using a PETSc application as an example. ExternalPetscSolverApp
 * create and destroys an external PETSc FEM/FDM solver.
 */
class ExternalPetscSolverApp : public MooseApp
{
public:
  ExternalPetscSolverApp(InputParameters parameters);
  virtual ~ExternalPetscSolverApp();

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);

#if LIBMESH_HAVE_PETSC
  /**
   * Return a time-stepping (TS) component that holds all the ingredients of applicaiton
   */
  TS & getExternalPETScTS() { return _ts; }
#endif
private:
#if LIBMESH_HAVE_PETSC
  /// Time-stepping (TS) object
  TS _ts;
#endif
};

