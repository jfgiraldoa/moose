//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// Moose Includes
#include "ThreadedElementLoop.h"

#include "libmesh/stored_range.h"

// Forward declarations
class DiracKernel;
template <typename T>
class MooseObjectTagWarehouse;
template <typename T>
class MooseObjectWarehouse;
class NonlinearSystemBase;

typedef StoredRange<std::set<const Elem *>::const_iterator, const Elem *> DistElemRange;

class ComputeDiracThread : public ThreadedElementLoop<DistElemRange>
{
public:
  ComputeDiracThread(FEProblemBase & feproblem, const std::set<TagID> & tags, bool _is_jacobian);

  // Splitting Constructor
  ComputeDiracThread(ComputeDiracThread & x, Threads::split);

  virtual ~ComputeDiracThread();

  virtual void subdomainChanged() override;
  virtual void pre() override;
  virtual void onElement(const Elem * elem) override;
  virtual void postElement(const Elem * /*elem*/) override;
  virtual void post() override;

  void join(const ComputeDiracThread & /*y*/);

protected:
  bool _is_jacobian;
  NonlinearSystemBase & _nl;

  const std::set<TagID> & _tags;

  /// Storage for DiracKernel objects
  MooseObjectTagWarehouse<DiracKernel> & _dirac_kernels;

  MooseObjectWarehouse<DiracKernel> * _dirac_warehouse;
};

