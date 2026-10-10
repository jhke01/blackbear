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

#include "GenericArrayNodalKernel.h"

/**
 * Dislocation-sink loss on one defect species' monomer component, using a dislocation density
 * that evolves in time via irradiation-induced climb (Pokor et al. 2004, Eqs 7-8) instead of the
 * fixed value the core void/loop kernels' own "sink" parameter assumes.
 *
 * rho(t) = (rho0^(-1/2) + K*b^2*t/2)^(-2) is a closed-form function of the current simulation
 * time alone -- it needs no coupling to any other variable and no tracked state beyond rho0, K,
 * and the dislocation network's own Burgers vector b. This kernel therefore replaces the role of
 * the core kernels' constant "sink" parameter (sink = D*Z*dislocation_density) with the
 * time-dependent D*Z*rho(t): one instance per species, each using that species' own D0/Em/Z but
 * the SAME rho0/K/b, matching the single, un-subscripted dislocation density Pokor's own Eqs
 * 3e/7/8 use for both species' sink terms.
 *
 * Unlike PointDefectRecombinationNodalKernel (which must couple to the OTHER species' array and
 * is therefore non-AD only), this kernel has no coupled variable at all, so it is AD-safe and
 * templated on is_ad like the core kernels.
 */
template <bool is_ad>
class GenericDislocationRecoveryNodalKernelTempl : public GenericArrayNodalKernel<is_ad>
{
public:
  static InputParameters validParams();
  GenericDislocationRecoveryNodalKernelTempl(const InputParameters & parameters);

protected:
  virtual void computeQpResidual(GenericRealEigenVector<is_ad> & residual) override;
  virtual void computeQpJacobian() override;

  /// Dislocation density at simulation time t, from the closed-form solution of Eq 8.
  Real rho(Real t) const;

  const Real _D0;
  const Real _Q;
  const Real _temperature;
  const Real _Z;
  const Real _rho0;
  const Real _inv_sqrt_rho0;
  const Real _K;
  const Real _burgers_vector;

  usingGenericArrayNodalKernelMembers;
};

typedef GenericDislocationRecoveryNodalKernelTempl<false> DislocationRecoveryNodalKernel;
typedef GenericDislocationRecoveryNodalKernelTempl<true> ADDislocationRecoveryNodalKernel;
