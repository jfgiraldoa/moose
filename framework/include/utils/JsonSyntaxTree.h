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
#include "FileLineInfo.h"
#include "json/json.h"
#include <string>
#include <vector>
#include <utility>

/**
 * Holds the syntax in a Json::Value tree
 */
class JsonSyntaxTree
{
public:
  JsonSyntaxTree(const std::string & search_string);
  virtual ~JsonSyntaxTree() {}

  /**
   * Add parameters to the tree
   * @param parent_path The parent syntax path that the action belongs to
   * @param path The path of the action
   * @param is_type Whether this belongs to a "<type>" or not
   * @param action Name of the action
   * @param is_action Wheter it is an action
   * @param params The InputParameters to add to the tree
   * @param lineinfo The FileLineInfo where the action/path was registered
   * @return Whether the parameters were added to the tree (ie if it matched the search string).
   */
  bool addParameters(const std::string & parent_path,
                     const std::string & path,
                     bool is_type,
                     const std::string & action,
                     bool is_action,
                     InputParameters * params,
                     const FileLineInfo & lineinfo,
                     const std::string & classname);

  /**
   * Add a task to the tree
   * @param path The path of the action
   * @param action Name of the action
   * @param task Name of the task
   * @param lineinfo The FileLineInfo where the action/task was registered
   */
  void addActionTask(const std::string & path,
                     const std::string & action,
                     const std::string & task,
                     const FileLineInfo & lineinfo);
  /**
   * Get the root of the tree.
   * @return The top level Json::Value holding the tree.
   */
  const moosecontrib::Json::Value & getRoot() const { return _root; }

  /**
   * Add an associated type to a block
   * @param path Path of the block
   * @param type Type name to associate the block with
   */
  void addSyntaxType(const std::string & path, const std::string type);

  /**
   * Add the global section to the output
   */
  void addGlobal();

protected:
  std::string buildOptions(const std::iterator_traits<InputParameters::iterator>::value_type & p,
                           bool & out_of_range_allowed);

  std::string prettyCppType(const std::string & cpp_type);
  std::string basicCppType(const std::string & cpp_type);
  size_t
  setParams(InputParameters * params, bool search_match, moosecontrib::Json::Value & all_params);

  std::string
  buildOutputString(const std::iterator_traits<InputParameters::iterator>::value_type & p);
  static std::vector<std::string> splitPath(const std::string & path);
  moosecontrib::Json::Value &
  getJson(const std::string & parent, const std::string & path, bool is_type);
  moosecontrib::Json::Value & getJson(const std::string & path);
  std::pair<std::string, std::string> getObjectLabel(const std::string & obj) const;
  std::pair<std::string, std::string> getActionLabel(const std::string & action) const;

  moosecontrib::Json::Value _root;
  std::string _search;

  ///@{
  /// Maps storing action/object name to the label and file location
  std::map<std::string, std::pair<std::string, std::string>> _action_label_map;
  std::map<std::string, std::pair<std::string, std::string>> _object_label_map;
  ///@}
};

