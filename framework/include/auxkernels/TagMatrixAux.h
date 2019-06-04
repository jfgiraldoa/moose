//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "TagVectorAux.h"

// Forward Declarations
class TagMatrixAux;

template <>
InputParameters validParams<TagMatrixAux>();

/**
 * For visualization or other purposes, the diagnal of the matrix of a tag
 * is extracted, and nodal values are assigned by using the matrix diagnal values.
 */
class TagMatrixAux : public AuxKernel
{
public:
  TagMatrixAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  TagID _tag_id;
  const VariableValue & _v;
};

