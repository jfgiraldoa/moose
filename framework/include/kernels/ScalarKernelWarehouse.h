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

#ifndef SCALARKERNELWAREHOUSE_H
#define SCALARKERNELWAREHOUSE_H

// MOOSE includes
#include "MooseObjectWarehouse.h"
#include "MooseTypes.h"
#include "MooseError.h"

// Forward declarations
class ScalarKernel;

/**
 * Holds kernels and provides some services
 */
class ScalarKernelWarehouse : public MooseObjectWarehouse<ScalarKernel>
{
public:
  ScalarKernelWarehouse(bool threaded = false);

  /**
   * Add Kernel to the storage structure
   */
  void addObject(std::shared_ptr<ScalarKernel> object, THREAD_ID tid = 0) override;

  ///@{
  /**
   * Methods for checking/getting variable kernels for a variable and SubdomainID
   */
  /*bool hasScalarVariableAllObjects(unsigned int variable_id,
                                     SubdomainID block_id,
                                     THREAD_ID tid = 0) const;*/
  const std::vector<std::shared_ptr<ScalarKernel>> & getScalarVariableAllObjects(
      unsigned int variable_id, /*SubdomainID block_id, */THREAD_ID tid = 0) const;
  ///@}

  /**
   * Update the active status of Kernels
   */
  //virtual void updateActive(THREAD_ID tid = 0) override;

protected:
  /// Variable based storage
  std::map<unsigned int, MooseObjectWarehouse<ScalarKernel>> _scalar_variable_kernel_storage;
};

#endif // SCALARKERNELWAREHOUSE_H
