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

#include "GenericDislocationRecoveryNodalKernel.h"

namespace
{
constexpr Real recovery_kB = 1.380649e-23;
constexpr Real recovery_eV = 1.602176634e-19;
} // namespace

registerMooseObject("BlackBearApp", DislocationRecoveryNodalKernel);
registerMooseObject("BlackBearApp", ADDislocationRecoveryNodalKernel);

template <bool is_ad>
InputParameters
GenericDislocationRecoveryNodalKernelTempl<is_ad>::validParams()
{
  InputParameters params = GenericArrayNodalKernel<is_ad>::validParams();
  params.addClassDescription(
      "Dislocation-sink loss on one defect species' monomer component (array component 0), "
      "using a dislocation density that evolves in time via irradiation-induced climb (Pokor et "
      "al. 2004, rho(t) = (rho0^(-1/2) + K*b^2*t/2)^(-2)) instead of a fixed value. Replaces the "
      "core void/loop kernels' own constant 'sink' parameter for species where this module is "
      "active; apply one instance per species, each with its own D0/Em_eV/Z but the same shared "
      "rho0/K/burgers_vector.");

  params.addRequiredParam<Real>("D0", "Diffusion prefactor D0 [m^2/s] for this species' monomer.");
  params.addRequiredParam<Real>("Em_eV", "Migration energy Em [eV] for this species' monomer.");
  params.addRequiredParam<Real>("temperature", "Temperature T [K].");
  params.addRequiredParam<Real>(
      "Z", "Dislocation-sink bias Z for this species' monomer (same quantity the core kernels' "
      "own constant-sink computation uses).");

  params.addParam<Real>(
      "rho0",
      "Initial dislocation density rho0 [m^-2] of the shared dislocation network. Required "
      "unless dislocation_density_postprocessor is given.");
  params.addParam<Real>(
      "K",
      "Recovery rate constant K [m^-1*s^-1] in d(rho)/dt = -K*b^2*rho^(3/2) (Pokor et al. 2004, "
      "Eq 8). Note: the paper's own text states K in bare m^-1, which is dimensionally "
      "inconsistent with this governing equation; treat K as a fit parameter in self-consistent "
      "SI units, not a literal copy of the paper's quoted numbers. Required unless "
      "dislocation_density_postprocessor is given.");
  params.addParam<PostprocessorName>(
      "dislocation_density_postprocessor",
      "Postprocessor that gives the dislocation density [m^-2] (for example "
      "CellularNetworkDislocationDensity). It replaces the closed form of Eq 8, so rho0 and K "
      "must not be given. The value is read as it stands at the start of each step.");
  params.addRequiredParam<Real>(
      "burgers_vector",
      "Burgers vector b [m] of the dislocation network itself -- distinct from any loop "
      "species' own geometric burgers_vector parameter.");

  return params;
}

template <bool is_ad>
GenericDislocationRecoveryNodalKernelTempl<is_ad>::GenericDislocationRecoveryNodalKernelTempl(
    const InputParameters & parameters)
  : GenericArrayNodalKernel<is_ad>(parameters),
    _D0(this->template getParam<Real>("D0")),
    _Q(this->template getParam<Real>("Em_eV") * recovery_eV),
    _temperature(this->template getParam<Real>("temperature")),
    _Z(this->template getParam<Real>("Z")),
    _rho0(this->isParamValid("rho0") ? this->template getParam<Real>("rho0") : 1.0),
    _inv_sqrt_rho0(1.0 / std::sqrt(_rho0)),
    _K(this->isParamValid("K") ? this->template getParam<Real>("K") : 0.0),
    _burgers_vector(this->template getParam<Real>("burgers_vector")),
    _rho_postprocessor(this->isParamValid("dislocation_density_postprocessor")
                           ? &this->getPostprocessorValue("dislocation_density_postprocessor")
                           : nullptr)
{
  if (_rho_postprocessor)
  {
    if (this->isParamValid("rho0") || this->isParamValid("K"))
      this->paramError("dislocation_density_postprocessor",
                       "cannot be combined with rho0 or K: the density comes from the postprocessor.");
  }
  else if (!this->isParamValid("rho0") || !this->isParamValid("K"))
    mooseError(this->type(), " requires rho0 and K unless dislocation_density_postprocessor is given.");

  if (_D0 <= 0.0)
    mooseError(this->type(), " requires D0 > 0.");
  if (_temperature <= 0.0)
    mooseError(this->type(), " requires temperature > 0.");
  if (_rho0 <= 0.0)
    mooseError(this->type(), " requires rho0 > 0.");
  if (_K < 0.0)
    mooseError(this->type(), " requires K >= 0.");
  if (_burgers_vector <= 0.0)
    mooseError(this->type(), " requires burgers_vector > 0.");
}

template <bool is_ad>
Real
GenericDislocationRecoveryNodalKernelTempl<is_ad>::rho(Real t) const
{
  if (_rho_postprocessor)
    return std::max(*_rho_postprocessor, 0.0);
  const Real inv_sqrt_rho =
      _inv_sqrt_rho0 + 0.5 * _K * _burgers_vector * _burgers_vector * t;
  return 1.0 / (inv_sqrt_rho * inv_sqrt_rho);
}

template <bool is_ad>
void
GenericDislocationRecoveryNodalKernelTempl<is_ad>::computeQpResidual(
    GenericRealEigenVector<is_ad> & residual)
{
  residual.resize(_count);
  residual.setZero();
  const GenericReal<is_ad> diffusivity = _D0 * std::exp(-_Q / (recovery_kB * _temperature));
  residual(0) = diffusivity * _Z * rho(this->_t) * _u[_qp](0);
}

template <>
void
GenericDislocationRecoveryNodalKernelTempl<false>::computeQpJacobian()
{
  const Real diffusivity = _D0 * std::exp(-_Q / (recovery_kB * _temperature));
  setJacobian(0, 0, diffusivity * _Z * rho(this->_t));
}

template <>
void
GenericDislocationRecoveryNodalKernelTempl<true>::computeQpJacobian()
{
  mooseError("Internal error: computeQpJacobian should never be called for the AD version");
}

template class GenericDislocationRecoveryNodalKernelTempl<false>;
template class GenericDislocationRecoveryNodalKernelTempl<true>;
