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

#include "GenericCascadeSpectrumNodalKernel.h"

registerMooseObject("BlackBearApp", CascadeSpectrumNodalKernel);
registerMooseObject("BlackBearApp", ADCascadeSpectrumNodalKernel);

template <bool is_ad>
InputParameters
GenericCascadeSpectrumNodalKernelTempl<is_ad>::validParams()
{
  InputParameters params = GenericArrayNodalKernel<is_ad>::validParams();
  params.addClassDescription(
      "Cascade-spectrum generation source on one defect species' explicit cluster sizes 2, 3, "
      "and 4 (array components 1, 2, 3), seeding clusters directly from the cascade instead of "
      "letting the entire generation rate go to the monomer (Pokor et al. 2004, Eqs 1a-1e). Each "
      "size n's source is generation*f_n/n, conserving total point-defect production regardless "
      "of f2/f3/f4.");

  params.addRequiredParam<Real>(
      "generation",
      "This species' shared monomer generation rate, generation = dpa_rate * "
      "cascade_efficiency (same quantity the core void/loop kernel's own 'generation' "
      "parameter uses).");
  params.addRequiredParam<Real>(
      "f2", "Fraction of generation seeded directly as cluster size 2.");
  params.addRequiredParam<Real>(
      "f3", "Fraction of generation seeded directly as cluster size 3.");
  params.addRequiredParam<Real>(
      "f4", "Fraction of generation seeded directly as cluster size 4.");

  return params;
}

template <bool is_ad>
GenericCascadeSpectrumNodalKernelTempl<is_ad>::GenericCascadeSpectrumNodalKernelTempl(
    const InputParameters & parameters)
  : GenericArrayNodalKernel<is_ad>(parameters),
    _generation(this->template getParam<Real>("generation")),
    _f2(this->template getParam<Real>("f2")),
    _f3(this->template getParam<Real>("f3")),
    _f4(this->template getParam<Real>("f4"))
{
  if (_f2 < 0.0 || _f3 < 0.0 || _f4 < 0.0)
    mooseError(this->type(), " requires f2, f3, and f4 to be nonnegative.");
  if (_f2 + _f3 + _f4 > 1.0)
    mooseError(this->type(), " requires f2 + f3 + f4 <= 1.");
  if (_count < 4)
    mooseError(this->type(),
               " requires the coupled array variable to represent at least cluster size 4, but "
               "it has only ",
               _count,
               " components.");
}

template <bool is_ad>
void
GenericCascadeSpectrumNodalKernelTempl<is_ad>::computeQpResidual(
    GenericRealEigenVector<is_ad> & residual)
{
  residual.resize(_count);
  residual.setZero();
  // Array component n - 1 holds explicit cluster size n (component 0 is the monomer). That holds
  // for sizes 2, 3, 4 only if the variable's explicit_max is at least 4; the kernel sees just the
  // component count, so the [DefectDynamics] action checks explicit_max and a hand-built input
  // has to make sure of it.
  residual(1) = -(_generation * _f2 / 2.0);
  residual(2) = -(_generation * _f3 / 3.0);
  residual(3) = -(_generation * _f4 / 4.0);
}

template <>
void
GenericCascadeSpectrumNodalKernelTempl<false>::computeQpJacobian()
{
  // The residual is a constant source with no dependence on the solution, so every Jacobian
  // entry is exactly zero -- nothing to set.
}

template <>
void
GenericCascadeSpectrumNodalKernelTempl<true>::computeQpJacobian()
{
  mooseError("Internal error: computeQpJacobian should never be called for the AD version");
}

template class GenericCascadeSpectrumNodalKernelTempl<false>;
template class GenericCascadeSpectrumNodalKernelTempl<true>;
