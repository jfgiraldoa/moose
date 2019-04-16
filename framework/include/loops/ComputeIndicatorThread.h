//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COMPUTEINDICATORTHREAD_H
#define COMPUTEINDICATORTHREAD_H

#include "ThreadedElementLoop.h"

#include "libmesh/elem_range.h"

// Forward declarations
class AuxiliarySystem;
class InternalSideIndicators;
class ExternalSideIndicators;

class ComputeIndicatorThread : public ThreadedElementLoop<ConstElemRange>
{
public:
  /**
   * @param fe_problem reference to the FEProblemBase we are computing on
   * @param sys reference to the AuxSystem we are computing on
   * @param indicator_whs Warehouse of Indicator objects.
   * @param finalize Whether or not we are just in the "finalize" stage or not.
   */
  ComputeIndicatorThread(FEProblemBase & fe_problem, bool finalize = false);

  // Splitting Constructor
  ComputeIndicatorThread(ComputeIndicatorThread & x, Threads::split split);

  virtual ~ComputeIndicatorThread();

  virtual void subdomainChanged() override;
  virtual void onElement(const Elem * elem) override;
  virtual void onBoundary(const Elem * elem, unsigned int side, BoundaryID bnd_id) override;
  virtual void onInternalSide(const Elem * elem, unsigned int side) override;
  virtual void postElement(const Elem * /*elem*/) override;
  virtual void post() override;

  void join(const ComputeIndicatorThread & /*y*/);

protected:
  FEProblemBase & _fe_problem;
  AuxiliarySystem & _aux_sys;

  /// Indicator Storage
  const MooseObjectWarehouse<Indicator> & _indicator_whs;

  /// InternalSideIndicator Storage
  const MooseObjectWarehouse<InternalSideIndicator> & _internal_side_indicators;

  /// ExternalSideIndicator Storage
  const MooseObjectWarehouse<ExternalSideIndicator> & _external_side_indicators;

  bool _finalize;
};

#endif // COMPUTEINDICATORTHREAD_H
