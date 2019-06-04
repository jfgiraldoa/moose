//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FEProblemBase.h"
#include "ParallelUniqueId.h"

template <typename RangeType, typename IteratorType>
class ThreadedNodeLoop
{
public:
  ThreadedNodeLoop(FEProblemBase & fe_problem);

  // Splitting Constructor
  ThreadedNodeLoop(ThreadedNodeLoop & x, Threads::split split);

  virtual ~ThreadedNodeLoop(){};

  void operator()(const RangeType & range);

  /**
   * Called before the node range loop
   */
  virtual void pre();

  /**
   * Called after the node range loop
   */
  virtual void post();

  /**
   * Called for each node
   */
  virtual void onNode(IteratorType & node_it);

  /**
   * Called after the node assembly is done (including surface assembling)
   *
   * @param node - active node
   */
  virtual void postNode(IteratorType & node_it);

  /**
   * Called if a MooseException is caught anywhere during the computation.
   * The single input parameter taken is a MooseException object.
   */
  virtual void caughtMooseException(MooseException & e)
  {
    std::string what(e.what());
    _fe_problem.setException(what);
  };

  /**
   * Whether or not the loop should continue.
   *
   * @return true to keep going, false to stop.
   */
  virtual bool keepGoing() { return !_fe_problem.hasException(); }

protected:
  FEProblemBase & _fe_problem;
  THREAD_ID _tid;
};

template <typename RangeType, typename IteratorType>
ThreadedNodeLoop<RangeType, IteratorType>::ThreadedNodeLoop(FEProblemBase & fe_problem)
  : _fe_problem(fe_problem)
{
}

template <typename RangeType, typename IteratorType>
ThreadedNodeLoop<RangeType, IteratorType>::ThreadedNodeLoop(ThreadedNodeLoop & x,
                                                            Threads::split /*split*/)
  : _fe_problem(x._fe_problem)
{
}

template <typename RangeType, typename IteratorType>
void
ThreadedNodeLoop<RangeType, IteratorType>::operator()(const RangeType & range)
{
  try
  {
    ParallelUniqueId puid;
    _tid = puid.id;

    pre();

    for (IteratorType nd = range.begin(); nd != range.end(); ++nd)
    {
      if (!keepGoing())
        break;

      onNode(nd);

      postNode(nd);
    }

    post();
  }
  catch (MooseException & e)
  {
    caughtMooseException(e);
  }
}

template <typename RangeType, typename IteratorType>
void
ThreadedNodeLoop<RangeType, IteratorType>::pre()
{
}

template <typename RangeType, typename IteratorType>
void
ThreadedNodeLoop<RangeType, IteratorType>::post()
{
}

template <typename RangeType, typename IteratorType>
void
ThreadedNodeLoop<RangeType, IteratorType>::onNode(IteratorType & /*node_it*/)
{
}

template <typename RangeType, typename IteratorType>
void
ThreadedNodeLoop<RangeType, IteratorType>::postNode(IteratorType & /*node_it*/)
{
}

