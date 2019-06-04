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
#include "MooseEnum.h"

// Forward declarations
class GeneratedMeshGenerator;

template <>
InputParameters validParams<GeneratedMeshGenerator>();

/**
 * Generates a line, square, or cube mesh with uniformly spaced or biased elements.
 */
class GeneratedMeshGenerator : public MeshGenerator
{
public:
  GeneratedMeshGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  /// The dimension of the mesh
  MooseEnum _dim;

  /// Number of elements in x, y, z direction
  unsigned int _nx, _ny, _nz;

  /// The min/max values for x,y,z component
  Real _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;

  /// All of the libmesh build_line/square/cube routines support an
  /// option to grade the mesh into the boundaries according to the
  /// spacing of the Gauss-Lobatto quadrature points.  Defaults to
  /// false, and cannot be used in conjunction with x, y, and z
  /// biasing.
  bool _gauss_lobatto_grid;

  /// The amount by which to bias the cells in the x,y,z directions.
  /// Must be in the range 0.5 <= _bias_x <= 2.0.
  /// _bias_x < 1 implies cells are shrinking in the x-direction.
  /// _bias_x==1 implies no bias (original mesh unchanged).
  /// _bias_x > 1 implies cells are growing in the x-direction.
  Real _bias_x, _bias_y, _bias_z;
};

