//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ThreadedElementLoop.h"
#include "MooseObjectTagWarehouse.h"

#include "libmesh/elem_range.h"

// Forward declarations
class FEProblemBase;
class NonlinearSystemBase;
class IntegratedBCBase;
class DGKernelBase;
class InterfaceKernelBase;
class Kernel;

class ComputeJacobianThread : public ThreadedElementLoop<ConstElemRange>
{
public:
  ComputeJacobianThread(FEProblemBase & fe_problem, const std::set<TagID> & tags);

  // Splitting Constructor
  ComputeJacobianThread(ComputeJacobianThread & x, Threads::split split);

  virtual ~ComputeJacobianThread();

  virtual void subdomainChanged() override;
  virtual void onElement(const Elem * elem) override;
  virtual void onBoundary(const Elem * elem, unsigned int side, BoundaryID bnd_id) override;
  virtual void onInternalSide(const Elem * elem, unsigned int side) override;
  virtual void onInterface(const Elem * elem, unsigned int side, BoundaryID bnd_id) override;
  virtual void postElement(const Elem * /*elem*/) override;
  virtual void post() override;

  void join(const ComputeJacobianThread & /*y*/);

protected:
  NonlinearSystemBase & _nl;

  unsigned int _num_cached;

  // Reference to BC storage structures
  MooseObjectTagWarehouse<IntegratedBCBase> & _integrated_bcs;

  MooseObjectWarehouse<IntegratedBCBase> * _ibc_warehouse;

  // Reference to DGKernel storage structure
  MooseObjectTagWarehouse<DGKernelBase> & _dg_kernels;

  MooseObjectWarehouse<DGKernelBase> * _dg_warehouse;

  // Reference to interface kernel storage structure
  MooseObjectTagWarehouse<InterfaceKernelBase> & _interface_kernels;

  MooseObjectWarehouse<InterfaceKernelBase> * _ik_warehouse;

  // Reference to Kernel storage structure
  MooseObjectTagWarehouse<KernelBase> & _kernels;

  // A pointer to different warehouse
  MooseObjectWarehouse<KernelBase> * _warehouse;

  /// Reference to ADKernel<JACOBIAN> storage structure
  MooseObjectTagWarehouse<KernelBase> & _ad_jacobian_kernels;

  /// Pointer to tags
  MooseObjectWarehouse<KernelBase> * _adjk_warehouse;

  const std::set<TagID> & _tags;

  virtual void computeJacobian();
  virtual void computeFaceJacobian(BoundaryID bnd_id);
  virtual void computeInternalFaceJacobian(const Elem * neighbor);
  virtual void computeInternalInterFaceJacobian(BoundaryID bnd_id);
};
