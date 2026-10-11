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
 * Cascade-spectrum generation source on one defect species' explicit cluster sizes 2, 3, and 4
 * (array components 1, 2, 3 -- component 0 is the monomer), seeding clusters directly from the
 * cascade instead of letting the entire generation rate go to the monomer (Pokor et al. 2004,
 * Eqs 1a-1e).
 *
 * To keep the total point-defect production rate conserved at exactly generation = dpa_rate *
 * cascade_efficiency regardless of how f2/f3/f4 split it (Pokor's own G(n) = generation*f_n form
 * manufactures n times too many point defects for cluster size n), each cluster size's source is
 * generation*f_n/n: a size-n cluster is injected at whatever rate keeps the point-defect count
 * conserved. This is a pure, state-independent constant, so the Jacobian contribution is exactly
 * zero and this kernel is AD-safe and templated on is_ad like the core kernels.
 *
 * This kernel does not reduce the monomer's own generation term; the enclosing
 * DefectDynamicsAction does that directly on the core void/loop kernel's own "generation"
 * parameter.
 */
template <bool is_ad>
class GenericCascadeSpectrumNodalKernelTempl : public GenericArrayNodalKernel<is_ad>
{
public:
  static InputParameters validParams();
  GenericCascadeSpectrumNodalKernelTempl(const InputParameters & parameters);

protected:
  virtual void computeQpResidual(GenericRealEigenVector<is_ad> & residual) override;
  virtual void computeQpJacobian() override;

  const Real _generation;
  const Real _f2;
  const Real _f3;
  const Real _f4;

  usingGenericArrayNodalKernelMembers;
};

typedef GenericCascadeSpectrumNodalKernelTempl<false> CascadeSpectrumNodalKernel;
typedef GenericCascadeSpectrumNodalKernelTempl<true> ADCascadeSpectrumNodalKernel;
