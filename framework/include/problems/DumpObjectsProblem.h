//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FEProblemBase.h"

class DumpObjectsProblem;
class DumpObjectsNonlinearSystem;

template <>
InputParameters validParams<FEProblem>();

/**
 * Specialization of SubProblem for dumping generated objects as input file syntax
 */
class DumpObjectsProblem : public FEProblemBase
{
public:
  DumpObjectsProblem(const InputParameters & parameters);

  void addVariable(const std::string & var_name,
                   const FEType & type,
                   Real scale_factor,
                   const std::set<SubdomainID> * const active_subdomains = NULL) override;
  void addScalarVariable(const std::string & var_name,
                         Order order,
                         Real scale_factor = 1.,
                         const std::set<SubdomainID> * const active_subdomains = NULL) override;
  void addAuxVariable(const std::string & var_name,
                      const FEType & type,
                      const std::set<SubdomainID> * const active_subdomains = NULL) override;
  void addAuxScalarVariable(const std::string & var_name,
                            Order order,
                            Real scale_factor = 1.,
                            const std::set<SubdomainID> * const active_subdomains = NULL) override;

  void addFunction(std::string type, const std::string & name, InputParameters parameters) override;

  void addKernel(const std::string & type,
                 const std::string & name,
                 InputParameters parameters) override;
  void addNodalKernel(const std::string & type,
                      const std::string & name,
                      InputParameters parameters) override;
  void addScalarKernel(const std::string & type,
                       const std::string & name,
                       InputParameters parameters) override;
  void addBoundaryCondition(const std::string & type,
                            const std::string & name,
                            InputParameters parameters) override;
  void addConstraint(const std::string & type,
                     const std::string & name,
                     InputParameters parameters) override;
  void addAuxKernel(const std::string & type,
                    const std::string & name,
                    InputParameters parameters) override;
  void addAuxScalarKernel(const std::string & type,
                          const std::string & name,
                          InputParameters parameters) override;
  void addDiracKernel(const std::string & type,
                      const std::string & name,
                      InputParameters parameters) override;
  void addDGKernel(const std::string & type,
                   const std::string & name,
                   InputParameters parameters) override;
  void addInterfaceKernel(const std::string & type,
                          const std::string & name,
                          InputParameters parameters) override;
  void addInitialCondition(const std::string & type,
                           const std::string & name,
                           InputParameters parameters) override;
  void addMaterial(const std::string & type,
                   const std::string & name,
                   InputParameters parameters) override;

  /// output input blocks for a given action path
  void dumpGeneratedSyntax(const std::string path);

  /// output data in solve
  virtual void solve() override;

  virtual void initialSetup() override {}
  virtual void advanceState() override {}
  virtual void timestepSetup() override {}
  virtual void execute(const ExecFlagType & /*exec_type*/) override {}
  virtual void outputStep(ExecFlagType /*type*/) override {}
  virtual void updateActiveObjects() override {}
  virtual void onTimestepEnd() override {}
  virtual void computeIndicators() override {}
  virtual void computeMarkers() override {}
  virtual bool adaptMesh() override { return false; }
  virtual void addLineSearch(const InputParameters & /*parameters*/) override {}

protected:
  void dumpObjectHelper(const std::string & system,
                        const std::string & type,
                        const std::string & name,
                        const InputParameters & parameters);

  void dumpVariableHelper(const std::string & system,
                          const std::string & var_name,
                          FEFamily family,
                          Order order,
                          Real scale_factor,
                          const std::set<SubdomainID> * const active_subdomains);

  /// build a text snippet of the minimal set of parameters that need to be specified
  std::string deduceNecessaryParameters(const std::string & type,
                                        const InputParameters & parameters);

  /// create a string map form parameter names to stringified parameter values
  std::map<std::string, std::string> stringifyParameters(const InputParameters & parameters);

  /// store input syntax to build objects generated by a specific action
  std::map<std::string, std::map<std::string, std::string>> _generated_syntax;

  std::shared_ptr<DumpObjectsNonlinearSystem> _nl_sys;
};

