//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_PETSC

#include "libmesh/petsc_nonlinear_solver.h"
#include <petscsnes.h>

class FEProblemBase;

/**
 * Wrapper of the libmesh ComputeLineSearchObject
 */
class ComputeLineSearchObjectWrapper
  : public libMesh::PetscNonlinearSolver<libMesh::Real>::ComputeLineSearchObject
{
public:
  ComputeLineSearchObjectWrapper(FEProblemBase & fe_problem);

  /**
   * Shim that calls into the MOOSE line search system using FEProblemBase::linesearch
   */
  void linesearch(SNESLineSearch line_search_object) override;

protected:
  FEProblemBase & _fe_problem;
};

#endif
