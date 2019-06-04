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

class GlobalParamsAction;

template <>
InputParameters validParams<GlobalParamsAction>();

class GlobalParamsAction : public Action
{
public:
  GlobalParamsAction(InputParameters params);

  virtual void act() override;

  /**
   * This function is here to remove parameters of a type so that global parameters
   * can potentially use the same named variable as a different type depending on the
   * application.
   */
  void remove(const std::string & name);

  template <typename T>
  inline T & setScalarParam(const std::string & name)
  {
    return parameters().set<T>(name);
  }

  template <typename T>
  inline std::vector<T> & setVectorParam(const std::string & name)
  {
    return parameters().set<std::vector<T>>(name);
  }

  template <typename T>
  inline std::vector<std::vector<T>> & setDoubleIndexParam(const std::string & name)
  {
    return parameters().set<std::vector<std::vector<T>>>(name);
  }
};
