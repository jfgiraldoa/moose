//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADMaterial.h"

#define usingStressUpdateBaseMembers                                                               \
  usingMaterialMembers;                                                                            \
  using ADStressUpdateBase<compute_stage>::updateState;                                            \
  using ADStressUpdateBase<compute_stage>::setQp;                                                  \
  using ADStressUpdateBase<compute_stage>::propagateQpStatefulProperties;                          \
  using ADStressUpdateBase<compute_stage>::requiresIsotropicTensor;                                \
  using ADStressUpdateBase<compute_stage>::computeTimeStepLimit;                                   \
  using ADStressUpdateBase<compute_stage>::_base_name

// Forward declarations
template <ComputeStage>
class ADStressUpdateBase;
template <typename>
class RankTwoTensorTempl;
typedef RankTwoTensorTempl<Real> RankTwoTensor;
typedef RankTwoTensorTempl<DualReal> DualRankTwoTensor;
template <typename>
class RankFourTensorTempl;
typedef RankFourTensorTempl<Real> RankFourTensor;
typedef RankFourTensorTempl<DualReal> DualRankFourTensor;
class InputParameters;

declareADValidParams(ADStressUpdateBase);

/**
 * ADStressUpdateBase is a material that is not called by MOOSE because
 * of the compute=false flag set in the parameter list.  This class is a base class
 * for materials that perform some internal computational
 * procedure (such as an iterative return-mapping procedure) to compute an
 * admissible state (which is usually an admissible stress that lies on or
 * within the yield surface, as well as a set of internal parameters such as
 * plastic strains).  The computational procedure must return the admissible stress
 * and a decomposition of the applied strain into elastic and inelastic components.
 * All materials inheriting from this class must be called by a separate material,
 * such as ComputeMultipleInelasticStress
 */
template <ComputeStage compute_stage>
class ADStressUpdateBase : public ADMaterial<compute_stage>
{
public:
  ADStressUpdateBase(const InputParameters & parameters);

  /**
   * Given a strain increment that results in a trial stress, perform some
   * procedure (such as an iterative return-mapping process) to produce
   * an admissible stress, an elastic strain increment and an inelastic
   * strain increment.
   * If _fe_problem.currentlyComputingJacobian() = true, then updateState also computes
   * d(stress)/d(strain) (or some approximation to it).
   *
   * This method is called by ComputeMultipleInelasticStress.
   * This method is pure virtual: all inheriting classes must overwrite this method.
   *
   * @param strain_increment Upon input: the strain increment.  Upon output: the elastic strain
   * increment
   * @param inelastic_strain_increment The inelastic_strain resulting from the iterative procedure
   * @param rotation_increment The finite-strain rotation increment
   * @param stress_new Upon input: the trial stress that results from applying strain_increment as
   * an elastic strain.  Upon output: the admissible stress
   * @param stress_old The old value of stress
   * @param elasticity_tensor The elasticity tensor
   */
  virtual void updateState(ADRankTwoTensor & strain_increment,
                           ADRankTwoTensor & inelastic_strain_increment,
                           const ADRankTwoTensor & rotation_increment,
                           ADRankTwoTensor & stress_new,
                           const RankTwoTensor & stress_old,
                           const ADRankFourTensor & elasticity_tensor,
                           const RankTwoTensor & elastic_strain_old) = 0;

  /// Sets the value of the global variable _qp for inheriting classes
  void setQp(unsigned int qp);

  /**
   * If updateState is not called during a timestep, this will be.  This method allows derived
   * classes to set internal parameters from their Old values, for instance
   */
  virtual void propagateQpStatefulProperties();

  /**
   * Does the model require the elasticity tensor to be isotropic?
   */
  virtual bool requiresIsotropicTensor() = 0;

  virtual Real computeTimeStepLimit();

  ///@{ Retained as empty methods to avoid a warning from Material.C in framework. These methods are unused in all inheriting classes and should not be overwritten.
  void resetQpProperties() final {}
  void resetProperties() final {}
  ///@}

protected:
  /// Name used as a prefix for all material properties related to the stress update model.
  const std::string _base_name;

  usingMaterialMembers;
};

