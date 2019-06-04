//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Action.h"
// libMesh
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"

class SetupQuadratureAction;

template <>
InputParameters validParams<SetupQuadratureAction>();

/**
 * Sets the quadrature
 */
class SetupQuadratureAction : public Action
{
public:
  SetupQuadratureAction(InputParameters parameters);

  virtual void act() override;

protected:
  QuadratureType _type;
  Order _order;
  Order _element_order;
  Order _side_order;
};

