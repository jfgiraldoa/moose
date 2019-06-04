//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MeshGenerator.h"

// Forward declarations
class FileMeshGenerator;

template <>
InputParameters validParams<FileMeshGenerator>();

/**
 * Generates a mesh by reading it from an file.
 */
class FileMeshGenerator : public MeshGenerator
{
public:
  FileMeshGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  const MeshFileName & _file_name;
};

