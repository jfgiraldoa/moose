//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralPostprocessor.h"

// Forward Declarations
class NumNodes;

template <>
InputParameters validParams<NumNodes>();

class NumNodes : public GeneralPostprocessor
{
public:
  NumNodes(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override {}

  /**
   * This will return the number of nodes in the system
   */
  virtual Real getValue() override;

private:
  const MeshBase & _mesh;
};

