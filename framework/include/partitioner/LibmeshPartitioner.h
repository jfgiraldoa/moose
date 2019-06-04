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
#include "MooseEnum.h"
#include "MoosePartitioner.h"

class LibmeshPartitioner;
class MooseMesh;

namespace libMesh
{
class SubdomainPartitioner;
}

template <>
InputParameters validParams<LibmeshPartitioner>();

class LibmeshPartitioner : public MoosePartitioner
{
public:
  LibmeshPartitioner(const InputParameters & params);
  virtual ~LibmeshPartitioner();

  virtual std::unique_ptr<Partitioner> clone() const;
  virtual void partition(MeshBase & mesh, const unsigned int n);
  virtual void partition(MeshBase & mesh);
  virtual void
  prepare_blocks_for_subdomain_partitioner(SubdomainPartitioner & subdomain_partitioner);

protected:
  virtual void _do_partition(MeshBase & mesh, const unsigned int n);

  std::unique_ptr<Partitioner> _partitioner;
  MooseEnum _partitioner_name;
  const std::vector<std::vector<SubdomainName>> & _subdomain_blocks;
  MooseMesh & _mesh;
};

