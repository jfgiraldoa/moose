//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InputParameters.h"
#include "ConsoleStreamInterface.h"
#include "Registry.h"
#include "PerfGraphInterface.h"

#include <string>
#include <ostream>

class Action;
class ActionWarehouse;
class ActionFactory;
class MooseMesh;
class FEProblemBase;
class Executioner;
class MooseApp;
class Factory;

template <>
InputParameters validParams<Action>();

/**
 * Base class for actions.
 */
class Action : public ConsoleStreamInterface, public PerfGraphInterface
{
public:
  Action(InputParameters parameters);

  virtual ~Action() {}

  /**
   * The method called externally that causes the action to act()
   */
  void timedAct();

protected:
  /**
   * Method to add a relationship manager for the objects being added to the system. Relationship
   * managers have to be added relatively early. In many cases before the Action::act() method
   * is called.
   * @param when_type The parameter indicating the normal time for adding either Geometric or
   *        Algebraic RelationshipManagers. It may not always be possible to add your
   *        RelationshipManager as early as you'd like. In these cases, your DistributedMesh may
   *        consume more memory during the problem setup.
   * @param moose_object_pars The MooseObject to inspect for RelationshipManagers to add
   */
  void addRelationshipManagers(Moose::RelationshipManagerType when_type,
                               const InputParameters & moose_object_pars);

public:
  /**
   * Method to add a relationship manager for the objects being added to the system. Relationship
   * managers have to be added relatively early. In many cases before the Action::act() method
   * is called.
   * @param when_type The parameter indicating the normal time for adding either Geometric or
   *        Algebraic RelationshipManagers. It may not always be possible to add your
   *        RelationshipManager as early as you'd like. In these cases, your DistributedMesh may
   *        consume more memory during the problem setup.
   */
  virtual void addRelationshipManagers(Moose::RelationshipManagerType when_type);

  /**
   * The name of the action
   */
  const std::string & name() const { return _name; }

  ///@{
  /**
   * Deprecated name methods, use name()
   */
  std::string getBaseName() const;
  std::string getShortName() const;
  ///@}

  const std::string & type() const { return _action_type; }

  InputParameters & parameters() { return _pars; }
  const InputParameters & parameters() const { return _pars; }

  const std::string & specificTaskName() const { return _specific_task_name; }

  const std::set<std::string> & getAllTasks() const { return _all_tasks; }

  ///@{
  /**
   * Retrieve a parameter for the object
   * @param name The name of the parameter
   * @return The value of the parameter
   */
  template <typename T>
  const T & getParam(const std::string & name) const;
  ///@}

  /**
   * Verifies that the requested parameter exists and is not NULL and returns it to the caller.
   * The template parameter must be a pointer or an error will be thrown.
   */
  template <typename T>
  T getCheckedPointerParam(const std::string & name, const std::string & error_string = "") const
  {
    return parameters().getCheckedPointerParam<T>(name, error_string);
  }

  inline bool isParamValid(const std::string & name) const { return _pars.isParamValid(name); }

  void appendTask(const std::string & task) { _all_tasks.insert(task); }

  /**
   * Emits an error prefixed with the file and line number of the given param (from the input
   * file) along with the full parameter path+name followed by the given args as the message.
   * If this object's parameters were not created directly by the Parser, then this function falls
   * back to the normal behavior of mooseError - only printing a message using the given args.
   */
  template <typename... Args>
  [[noreturn]] void paramError(const std::string & param, Args... args) {
    auto prefix = param + ": ";
    if (!_pars.inputLocation(param).empty())
      prefix = _pars.inputLocation(param) + ": (" + _pars.paramFullpath(param) + "):\n";
    mooseError(prefix, args...);
  }

  /**
   * Emits a warning prefixed with the file and line number of the given param (from the input
   * file) along with the full parameter path+name followed by the given args as the message.
   * If this object's parameters were not created directly by the Parser, then this function falls
   * back to the normal behavior of mooseWarning - only printing a message using the given args.
   */
  template <typename... Args>
  void paramWarning(const std::string & param, Args... args)
  {
    auto prefix = param + ": ";
    if (!_pars.inputLocation(param).empty())
      prefix = _pars.inputLocation(param) + ": (" + _pars.paramFullpath(param) + "):\n";
    mooseWarning(prefix, args...);
  }

  /**
   * Emits an informational message prefixed with the file and line number of the given param
   * (from the input file) along with the full parameter path+name followed by the given args as
   * the message.  If this object's parameters were not created directly by the Parser, then this
   * function falls back to the normal behavior of mooseInfo - only printing a message using
   * the given args.
   */
  template <typename... Args>
  void paramInfo(const std::string & param, Args... args)
  {
    auto prefix = param + ": ";
    if (!_pars.inputLocation(param).empty())
      prefix = _pars.inputLocation(param) + ": (" + _pars.paramFullpath(param) + "):\n";
    mooseInfo(prefix, args...);
  }

protected:
  /**
   * Method to add objects to the simulation or perform other setup tasks.
   */
  virtual void act() = 0;

  /// Input parameters for the action
  InputParameters _pars;

  // The registered syntax for this block if any
  std::string _registered_identifier;

  /// The name of the action
  std::string _name;

  // The type name of this Action instance
  std::string _action_type;

  /// The MOOSE application this is associated with
  MooseApp & _app;

  /// The Factory associated with the MooseApp
  Factory & _factory;

  /// Builds Actions
  ActionFactory & _action_factory;

  /**
   * This member will only be populated if this Action instance is only designed to
   * handle one task.  This happens when an Action is registered with several pieces
   * of syntax in which case separate instances are built to handle the different
   * incoming parameter values.
   */
  std::string _specific_task_name;

  /**
   * A list of all the tasks that this Action will satisfy.
   * Note: That this is _not_ populated at construction time.  However, all tasks will be
   *       added prior to act().
   */
  std::set<std::string> _all_tasks;

  /// Reference to ActionWarehouse where we store object build by actions
  ActionWarehouse & _awh;

  /// The current action (even though we have seperate instances for each action)
  const std::string & _current_task;

  std::shared_ptr<MooseMesh> & _mesh;
  std::shared_ptr<MooseMesh> & _displaced_mesh;

  /// Convenience reference to a problem this action works on
  std::shared_ptr<FEProblemBase> & _problem;

  /// Timers
  PerfID _act_timer;
};

template <typename T>
const T &
Action::getParam(const std::string & name) const
{
  return InputParameters::getParamHelper(name, _pars, static_cast<T *>(0));
}

