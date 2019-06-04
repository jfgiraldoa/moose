//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseObjectUnitTest.h"
#include "HeliumFluidProperties.h"

class HeliumFluidPropertiesTest : public MooseObjectUnitTest
{
public:
  HeliumFluidPropertiesTest() : MooseObjectUnitTest("MooseUnitApp") { buildObjects(); }

protected:
  void buildObjects()
  {
    InputParameters uo_pars = _factory.getValidParams("HeliumFluidProperties");
    _fe_problem->addUserObject("HeliumFluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObject<HeliumFluidProperties>("fp");
  }

  const HeliumFluidProperties * _fp;
};
