//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseMesh.h"
#include "ThreadedNodeLoop.h"
#include "MooseObjectTagWarehouse.h"

class AuxiliarySystem;
class NodalKernel;

class ComputeNodalKernelBCJacobiansThread
  : public ThreadedNodeLoop<ConstBndNodeRange, ConstBndNodeRange::const_iterator>
{
public:
  ComputeNodalKernelBCJacobiansThread(FEProblemBase & fe_problem,
                                      MooseObjectTagWarehouse<NodalKernel> & nodal_kernels,
                                      const std::set<TagID> & tags);

  // Splitting Constructor
  ComputeNodalKernelBCJacobiansThread(ComputeNodalKernelBCJacobiansThread & x,
                                      Threads::split split);

  virtual void pre() override;

  virtual void onNode(ConstBndNodeRange::const_iterator & node_it) override;

  void join(const ComputeNodalKernelBCJacobiansThread & /*y*/);

protected:
  FEProblemBase & _fe_problem;

  AuxiliarySystem & _aux_sys;

  const std::set<TagID> & _tags;

  MooseObjectTagWarehouse<NodalKernel> & _nodal_kernels;

  MooseObjectWarehouse<NodalKernel> * _nkernel_warehouse;

  /// Number of contributions cached up
  unsigned int _num_cached;
};

