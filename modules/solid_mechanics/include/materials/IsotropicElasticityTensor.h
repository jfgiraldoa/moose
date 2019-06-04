//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElasticityTensor.h"

/**
 * Defines an Isotropic Elasticity Tensor.
 *
 * The input is any two of the following:
 *
 * Youngs Modulus
 * Poissons Ration
 * First and Second Lame Coefficients (lambda and mu)
 * Bulk Modulus
 *
 * Internally this class uses the Lame Coefficients.
 * Everything is is transformed to these coefficients.
 *
 * Note that by default this tensor is _constant_... meaning
 * that it never changes after the first time it is computed.
 *
 * If you want to modify this behavior you can pass in
 * false to the constructor.
 */
class IsotropicElasticityTensor : public ElasticityTensor
{
public:
  IsotropicElasticityTensor(const bool constant = true);

  /**
   * Set the first Lame Coefficient.
   */
  void setLambda(const Real lambda);

  /**
   * Set the second Lame Coefficient.
   */
  void setMu(const Real mu);

  /**
   * Set the Young's Modulus
   */
  void setYoungsModulus(const Real E);

  /**
   * Set Poissons Ratio
   */
  void setPoissonsRatio(const Real nu);

  /**
   * Set the Bulk Modulus
   */
  void setBulkModulus(const Real k);

  /**
   * Set the shear modulus... same thing as Mu
   */
  void setShearModulus(const Real k);

  virtual ~IsotropicElasticityTensor() {}

protected:
  bool _lambda_set, _mu_set, _E_set, _nu_set, _k_set;

  Real _lambda, _mu, _E, _nu, _k;

  /**
   * Fill in the matrix.
   */
  virtual void calculateEntries(unsigned int qp);

  /**
   * Calculates lambda and mu based on what has been set.
   *
   * These are based on Michael Tonks's's notes
   */
  void calculateLameCoefficients();

  /**
   * Computes a single entry of C_ijkl.
   *
   * Note that the formula for this came from Michael Tonks on page 234 of his notes.
   */
  Real isotropicEntry(const unsigned int i, const unsigned j, const unsigned k, const unsigned l);
};

