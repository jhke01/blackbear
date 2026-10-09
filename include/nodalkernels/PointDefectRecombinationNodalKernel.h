/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*                       BlackBear                              */
/*                                                              */
/*           (c) 2017 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#pragma once

#include "ArrayNodalKernel.h"

/**
 * Point-defect recombination loss on one defect family's monomer component,
 * coupled to the other family's own array variable.
 *
 * Recombination consumes one vacancy monomer and one interstitial monomer
 * per event, at rate recombination_rate * C_v1 * C_i1
 * (recombination_rate = 4*pi*r_rec*(D_i+D_v)/V_at, matching
 * DefectDynamicsAMSSInput/DefectDynamics_Sundials.cpp). This kernel adds that
 * loss to the monomer component (component 0) of the variable it is applied
 * to; two instances, one on each family's array variable with "coupled_var"
 * set to the other, are needed to represent recombination symmetrically.
 *
 * This kernel is intentionally non-AD only. MOOSE's only array-variable
 * coupling accessor, Coupleable::coupledArrayValue(), returns a plain
 * ArrayVariableValue with no AD derivative tracking, and
 * ADArrayNodalKernel::setJacobian() unconditionally errors, so an AD variant
 * could never actually assemble the cross-variable Jacobian block: it would
 * silently leave it zero, degrading Newton's convergence rate without any
 * compile- or run-time indication of the problem. The residual (and
 * therefore the converged solution) does not depend on this limitation; only
 * the ability to differentiate automatically with respect to the *coupled*
 * variable does.
 */
class PointDefectRecombinationNodalKernel : public ArrayNodalKernel
{
public:
  static InputParameters validParams();
  PointDefectRecombinationNodalKernel(const InputParameters & parameters);

protected:
  virtual void computeQpResidual(RealEigenVector & residual) override;
  virtual void computeQpJacobian() override;
  virtual void computeQpOffDiagJacobian(unsigned int jvar) override;

  const Real _recombination_rate;
  const ArrayVariableValue & _coupled_monomer;
  const unsigned int _coupled_var_num;
};
