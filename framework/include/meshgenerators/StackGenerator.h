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
#include "libmesh/replicated_mesh.h"
#include "MooseEnum.h"

// Forward declarations
class StackGenerator;

template <>
InputParameters validParams<StackGenerator>();

/**
 * Take several 3D meshes and stitch them on top of each other like a stack.
 */
class StackGenerator : public MeshGenerator
{
public:
  StackGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  /// The dimension of the mesh
  MooseEnum _dim;

  /// The meshgenerators to read
  const std::vector<MeshGeneratorName> & _input_names;

  /// Height (z) of the bottom of the final mesh
  const Real _bottom_height;

  // Holds pointers to the pointers to the meshes.
  std::vector<std::unique_ptr<MeshBase> *> _mesh_ptrs;

  /// The meshes to be stitched together.
  std::vector<std::unique_ptr<ReplicatedMesh>> _meshes;

  Real computeWidth(const MeshBase & mesh, const int & dim);
};

