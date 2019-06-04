//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "SystemBase.h"

// Forward declarations
class DisplacedProblem;
namespace libMesh
{
class ExplicitSystem;
template <typename>
class TransientSystem;
typedef TransientSystem<ExplicitSystem> TransientExplicitSystem;
}

class DisplacedSystem : public SystemBase
{
public:
  DisplacedSystem(DisplacedProblem & problem,
                  SystemBase & undisplaced_system,
                  const std::string & name,
                  Moose::VarKindType var_kind);
  virtual ~DisplacedSystem();

  virtual void init() override;

  virtual NumericVector<Number> & getVector(TagID tag_id) override
  {
    return _undisplaced_system.getVector(tag_id);
  }
  virtual const NumericVector<Number> & getVector(TagID tag_id) const override
  {
    return _undisplaced_system.getVector(tag_id);
  }

  virtual TagID residualVectorTag() override { return _undisplaced_system.residualVectorTag(); }

  virtual TagID systemMatrixTag() override { return _undisplaced_system.systemMatrixTag(); }

  virtual TagID timeVectorTag() override { return _undisplaced_system.timeVectorTag(); }

  virtual TagID nonTimeVectorTag() override { return _undisplaced_system.nonTimeVectorTag(); }

  virtual void associateVectorToTag(NumericVector<Number> & vec, TagID tag) override
  {
    _undisplaced_system.associateVectorToTag(vec, tag);
  }

  virtual void disassociateVectorFromTag(NumericVector<Number> & vec, TagID tag) override
  {
    _undisplaced_system.disassociateVectorFromTag(vec, tag);
  }

  virtual void disassociateAllTaggedVectors() override
  {
    _undisplaced_system.disassociateAllTaggedVectors();
  }

  virtual void associateMatrixToTag(SparseMatrix<Number> & matrix, TagID tag) override
  {
    _undisplaced_system.associateMatrixToTag(matrix, tag);
  }

  virtual void disassociateMatrixFromTag(SparseMatrix<Number> & matrix, TagID tag) override
  {
    _undisplaced_system.disassociateMatrixFromTag(matrix, tag);
  }

  virtual void disassociateAllTaggedMatrices() override
  {
    _undisplaced_system.disassociateAllTaggedMatrices();
  }

  virtual NumericVector<Number> & getVector(const std::string & name) override;

  virtual NumericVector<Number> & serializedSolution() override
  {
    return _undisplaced_system.serializedSolution();
  }

  const NumericVector<Number> * const & currentSolution() const override
  {
    return _undisplaced_system.currentSolution();
  }

  NumericVector<Number> & solution() override { return _undisplaced_system.solution(); }
  NumericVector<Number> & solutionOld() override;
  NumericVector<Number> & solutionOlder() override;
  NumericVector<Number> * solutionPreviousNewton() override { return NULL; }

  const NumericVector<Number> & solution() const override { return _undisplaced_system.solution(); }
  const NumericVector<Number> & solutionOld() const override;
  const NumericVector<Number> & solutionOlder() const override;
  const NumericVector<Number> * solutionPreviousNewton() const override { return NULL; }

  NumericVector<Number> * solutionUDot() override { return _undisplaced_system.solutionUDot(); }
  NumericVector<Number> * solutionUDotDot() override
  {
    return _undisplaced_system.solutionUDotDot();
  }
  NumericVector<Number> * solutionUDotOld() override
  {
    return _undisplaced_system.solutionUDotOld();
  }
  NumericVector<Number> * solutionUDotDotOld() override
  {
    return _undisplaced_system.solutionUDotDotOld();
  }
  const NumericVector<Number> * solutionUDot() const override
  {
    return _undisplaced_system.solutionUDot();
  }
  const NumericVector<Number> * solutionUDotDot() const override
  {
    return _undisplaced_system.solutionUDotDot();
  }
  const NumericVector<Number> * solutionUDotOld() const override
  {
    return _undisplaced_system.solutionUDotOld();
  }
  const NumericVector<Number> * solutionUDotDotOld() const override
  {
    return _undisplaced_system.solutionUDotDotOld();
  }

  virtual Number & duDotDu() override { return _undisplaced_system.duDotDu(); }
  virtual Number & duDotDotDu() override { return _undisplaced_system.duDotDotDu(); }
  virtual const Number & duDotDu() const override { return _undisplaced_system.duDotDu(); }
  virtual const Number & duDotDotDu() const override { return _undisplaced_system.duDotDotDu(); }

  /**
   * Return the residual copy from the NonlinearSystem
   * @return Residual copy
   */
  virtual NumericVector<Number> & residualCopy() override
  {
    return _undisplaced_system.residualCopy();
  }
  virtual NumericVector<Number> & residualGhosted() override
  {
    return _undisplaced_system.residualGhosted();
  }

  virtual void augmentSendList(std::vector<dof_id_type> & send_list) override
  {
    _undisplaced_system.augmentSendList(send_list);
  }

  /**
   * This is an empty function since the displaced system doesn't have a matrix!
   * All sparsity pattern modification will be taken care of by the undisplaced system directly
   */
  virtual void augmentSparsity(SparsityPattern::Graph & /*sparsity*/,
                               std::vector<dof_id_type> & /*n_nz*/,
                               std::vector<dof_id_type> & /*n_oz*/) override
  {
  }

  /**
   * Adds this variable to the list of variables to be zeroed during each residual evaluation.
   * @param var_name The name of the variable to be zeroed.
   */
  virtual void addVariableToZeroOnResidual(std::string var_name) override
  {
    _undisplaced_system.addVariableToZeroOnResidual(var_name);
  }

  /**
   * Adds this variable to the list of variables to be zeroed during each jacobian evaluation.
   * @param var_name The name of the variable to be zeroed.
   */
  virtual void addVariableToZeroOnJacobian(std::string var_name) override
  {
    _undisplaced_system.addVariableToZeroOnJacobian(var_name);
  }

  /**
   * Zero out the solution for the list of variables passed in.
   */
  virtual void zeroVariables(std::vector<std::string> & vars_to_be_zeroed) override
  {
    _undisplaced_system.zeroVariables(vars_to_be_zeroed);
  }

  virtual bool hasVector(TagID tag_id) const override
  {
    return _undisplaced_system.hasVector(tag_id);
  }

  virtual bool hasMatrix(TagID tag_id) const override
  {
    return _undisplaced_system.hasMatrix(tag_id);
  }

  virtual SparseMatrix<Number> & getMatrix(TagID tag) override
  {
    return _undisplaced_system.getMatrix(tag);
  }
  virtual const SparseMatrix<Number> & getMatrix(TagID tag) const override
  {
    return _undisplaced_system.getMatrix(tag);
  }

  virtual TransientExplicitSystem & sys() { return _sys; }

  virtual System & system() override;
  virtual const System & system() const override;

  using SystemBase::addTimeIntegrator;
  void addTimeIntegrator(std::shared_ptr<TimeIntegrator> ti) override;

protected:
  SystemBase & _undisplaced_system;
  TransientExplicitSystem & _sys;
};

