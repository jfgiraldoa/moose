/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "ScalarKernelWarehouse.h"

// MOOSE includes
#include "ScalarKernel.h"
#include "MooseVariableScalar.h"
#include "TimeDerivative.h"
#include "TimeKernel.h"

ScalarKernelWarehouse::ScalarKernelWarehouse(bool threaded) : MooseObjectWarehouse<ScalarKernel>(threaded) {}

void
ScalarKernelWarehouse::addObject(std::shared_ptr<ScalarKernel> object, THREAD_ID tid)
{
  // Add object to the general storage
  MooseObjectWarehouse<ScalarKernel>::addObject(object, tid);

  // Add object to the variable based storage
  _scalar_variable_kernel_storage[object->variable().number()].addObject(object, tid);
}
/*
bool
ScalarKernelWarehouse::hasScalarVariableAllObjects(unsigned int variable_id,
                                               //SubdomainID block_id,
                                               THREAD_ID tid) const
{
  checkThreadID(tid);
  std::map<unsigned int, MooseObjectWarehouse<ScalarKernel>>::const_iterator iter =
      _scalar_variable_kernel_storage.find(variable_id);
  return (iter != _scalar_variable_kernel_storage.end() &&
          //iter->second.hasActiveBlockObjects(block_id, tid));
          iter->second.hasObjects(tid));
}*/

const std::vector<std::shared_ptr<ScalarKernel>> &
ScalarKernelWarehouse::getScalarVariableAllObjects(unsigned int variable_id,
                                               //SubdomainID block_id,
                                               THREAD_ID tid) const
{
  checkThreadID(tid);
  const auto iter = _scalar_variable_kernel_storage.find(variable_id);
  mooseAssert(iter != _scalar_variable_kernel_storage.end(),
              "Unable to located variable kernels for the given variable id: " << variable_id
                                                                               << ".");
  return iter->second.getObjects(tid);
  //return iter->second.getActiveBlockObjects(block_id, tid);
}

/*
void
ScalarKernelWarehouse::updateActive(THREAD_ID tid)
{
  MooseObjectWarehouse<ScalarKernel>::updateActive(tid);

  for (auto & it : _scalar_variable_kernel_storage)
    it.second.updateActive(tid);
}
*/
