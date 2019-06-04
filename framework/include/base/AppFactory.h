//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include <vector>

#include "MooseApp.h"

// Forward declarations
class InputParameters;

/**
 * Macros
 */
#define registerApp(name) AppFactory::instance().reg<name>(#name)

/**
 * alias to wrap shared pointer type
 */
using MooseAppPtr = std::shared_ptr<MooseApp>;

/**
 * alias for validParams function
 */
using paramsPtr = InputParameters (*)();

/**
 * alias for method to build objects
 */
using appBuildPtr = MooseAppPtr (*)(const InputParameters & parameters);

/**
 * alias for registered Object iterator
 */
using registeredMooseAppIterator = std::map<std::string, paramsPtr>::iterator;

/**
 * Build an object of type T
 */
template <class T>
MooseAppPtr
buildApp(const InputParameters & parameters)
{
  return std::make_shared<T>(parameters);
}

/**
 * Generic AppFactory class for building Application objects
 */
class AppFactory
{
public:
  /**
   * Get the instance of the AppFactory
   * @return Pointer to the AppFactory instance
   */
  static AppFactory & instance();

  virtual ~AppFactory();

  /**
   * Helper function for creating a MooseApp from command-line arguments.
   */
  static MooseAppPtr createAppShared(const std::string & app_type,
                                     int argc,
                                     char ** argv,
                                     MPI_Comm comm_word = MPI_COMM_WORLD);

  /**
   * Register a new object
   * @param name Name of the object to register
   */
  template <typename T>
  void reg(const std::string & name)
  {
    if (_name_to_build_pointer.find(name) == _name_to_build_pointer.end())
    {
      _name_to_build_pointer[name] = &buildApp<T>;
      _name_to_params_pointer[name] = &validParams<T>;
    }
  }

  /**
   * Get valid parameters for the object
   * @param name Name of the object whose parameter we are requesting
   * @return Parameters of the object
   */
  InputParameters getValidParams(const std::string & name);

  /**
   * Build an application object (must be registered)
   * @param app_type Type of the application being constructed
   * @param name Name for the object
   * @param parameters Parameters this object should have
   * @return The created object
   */
  MooseAppPtr createShared(const std::string & app_type,
                           const std::string & name,
                           InputParameters parameters,
                           MPI_Comm COMM_WORLD_IN);

  ///@{
  /**
   * Returns iterators to the begin/end of the registered objects data structure: a name ->
   * validParams function pointer.
   */
  registeredMooseAppIterator registeredObjectsBegin() { return _name_to_params_pointer.begin(); }
  registeredMooseAppIterator registeredObjectsEnd() { return _name_to_params_pointer.end(); }
  ///@}

  /**
   * Returns a Boolean indicating whether an application type has been registered
   */
  bool isRegistered(const std::string & app_name) const;

protected:
  std::map<std::string, appBuildPtr> _name_to_build_pointer;

  std::map<std::string, paramsPtr> _name_to_params_pointer;

  static AppFactory _instance;

private:
  // Private constructor for singleton pattern
  AppFactory() {}
};

