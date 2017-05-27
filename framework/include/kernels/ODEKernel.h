/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef ODEKERNEL_H
#define ODEKERNEL_H

#include "ScalarKernel.h"

// Forward Declarations
class ODEKernel;

template <>
InputParameters validParams<ODEKernel>();

/**
 *
 */
class ODEKernel : public ScalarKernel
{
public:
  ODEKernel(const InputParameters & parameters);

  virtual void reinit() override;
  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  virtual void computeOffDiagJacobian(unsigned int jvar) override;

protected:
  virtual Real computeQpResidual() = 0;
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  //virtual void precalculateOffDiagJacobian(unsigned int jvar);
  Assembly & _assembly;
  unsigned int _qp;
  QBase *& _qrule;
  const MooseArray<Real> & _JxW;
  const MooseArray<Real> & _coord;
  const VariablePhiValue & _phi;
};

#endif /* ODEKERNEL_H */
