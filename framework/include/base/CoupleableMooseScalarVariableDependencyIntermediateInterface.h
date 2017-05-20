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

#ifndef COUPLEABLEMOOSESCALARVARIABLEDEPENDENCYINTERMEDIATEINTERFACE_H
#define COUPLEABLEMOOSESCALARVARIABLEDEPENDENCYINTERMEDIATEINTERFACE_H

#include "Coupleable.h"
#include "ScalarCoupleable.h"
//#include "MooseVariableInterface.h"
#include "MooseVariableDependencyInterface.h"

/**
 * Intermediate base class that ties together all the interfaces for getting
 * MooseVariables with the MooseVariableDependencyInterface
 */
class CoupleableMooseScalarVariableDependencyIntermediateInterface
    : public Coupleable,
      public ScalarCoupleable,
      //public MooseVariableInterface,
      public MooseVariableDependencyInterface
{
public:
  CoupleableMooseScalarVariableDependencyIntermediateInterface(const MooseObject * moose_object,
                                                         bool nodal)
    : Coupleable(moose_object, nodal),
      ScalarCoupleable(moose_object)
      //MooseVariableInterface(moose_object, nodal)
  {
    const std::vector<MooseVariable *> & coupled_vars = getCoupledMooseVars();
    for (unsigned int i = 0; i < coupled_vars.size(); i++)
      addMooseVariableDependency(coupled_vars[i]);

    //addMooseVariableDependency(mooseVariable());
  }
};

#endif // COUPLEABLEMOOSESCALARVARIABLEDEPENDENCYINTERMEDIATEINTERFACE_H
