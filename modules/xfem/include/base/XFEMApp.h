//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseApp.h"

class XFEMApp;

template <>
InputParameters validParams<XFEMApp>();

class XFEMApp : public MooseApp
{
public:
  XFEMApp(const InputParameters & parameters);
  virtual ~XFEMApp();

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
  static void registerObjectDepends(Factory & factory);
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
  static void registerExecFlags(Factory & factory);
  static void associateSyntaxDepends(Syntax & syntax, ActionFactory & action_factory);
};

