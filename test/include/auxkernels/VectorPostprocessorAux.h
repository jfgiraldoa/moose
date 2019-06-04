//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

// Forward Declarations
class VectorPostprocessorAux;

template <>
InputParameters validParams<VectorPostprocessorAux>();

/**
 * Coupled auxiliary value
 */
class VectorPostprocessorAux : public AuxKernel
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  VectorPostprocessorAux(const InputParameters & parameters);

  virtual ~VectorPostprocessorAux() {}

protected:
  virtual Real computeValue();

  const VectorPostprocessorValue & _vpp;

  unsigned int _index;
};

