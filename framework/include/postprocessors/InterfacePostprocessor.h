//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceUserObject.h"
#include "Postprocessor.h"

// Forward Declarations
class InterfacePostprocessor;

/**
 * BASE CLASS FOR IMPLEMENTING INTERFACE POSTRPOCESSOR
 * ALL INTERFACE POSTPROCESOR SHOULD DERIVE FROM THIS CLASS
 * THIS CALSS ALASO COMPUTE THE INTERFACE AREA
 */

template <>
InputParameters validParams<InterfacePostprocessor>();

class InterfacePostprocessor : public InterfaceUserObject, public Postprocessor
{
public:
  InterfacePostprocessor(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & y) override;
  /**
   * This is called _after_ execute() and _after_ threadJoin()!  This is probably where you want to
   * do MPI communication!
   * Finalize is not required for Postprocessor implementations since work may be done in
   * getValue().
   */
  virtual Real getValue() override
  {
    // compute the interface area based on the are of the master side
    // this come at no cost so why not mek it avaialable for all subclasses
    gatherSum(_interface_master_area);
    return 0;
  }

  // becasue getValue already does the work we just override to nothing
  virtual void finalize() override {}

protected:
  /// the area of the master side of the interface
  Real _interface_master_area;
};

