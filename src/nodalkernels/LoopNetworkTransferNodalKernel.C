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

#include "LoopNetworkTransferNodalKernel.h"

#include <cmath>

registerMooseObject("BlackBearApp", LoopNetworkTransferNodalKernel);

namespace
{
constexpr Real transfer_kB_eV = 8.617333262e-5;
}

InputParameters
LoopNetworkTransferNodalKernel::validParams()
{
  InputParameters params = ArrayNodalKernel::validParams();
  params.addClassDescription(
      "Loss of interstitial loops to the dislocation network at the rate k = v_l / L, with the "
      "climb velocity v_l = 2 Z_iv D_i c_i / b and the cell size L = C_A C_rho sqrt(pi / (rho f)) "
      "of Rest (2006), the counterpart of the line-length source of "
      "CellularNetworkDislocationDensity.");
  params.addRequiredParam<PostprocessorName>("dislocation_density_postprocessor",
                                             "Network density rho [m^-2].");
  params.addRequiredParam<Real>("D0", "Interstitial diffusion prefactor D0 [m^2/s].");
  params.addRequiredParam<Real>("Em_eV", "Interstitial migration energy [eV].");
  params.addRequiredParam<Real>("temperature", "Temperature [K].");
  params.addRequiredParam<Real>("burgers_vector", "Burgers vector b [m].");
  params.addRequiredParam<Real>("bias", "Relative interstitial/vacancy bias Z_iv.");
  params.addParam<Real>("ca_crho", 2.0, "Product C_A C_rho of the cell geometry.");
  params.addParam<Real>("poisson_ratio", 0.3, "Poisson ratio nu of f(nu) = (1 - nu/2)/(1 - nu).");
  params.addRangeCheckedParam<unsigned int>(
      "minimum_size", 2, "minimum_size >= 2", "Smallest loop size that can join the network.");
  return params;
}

LoopNetworkTransferNodalKernel::LoopNetworkTransferNodalKernel(const InputParameters & parameters)
  : ArrayNodalKernel(parameters),
    _coefficient(2.0 * getParam<Real>("bias") * getParam<Real>("D0") *
                 std::exp(-getParam<Real>("Em_eV") / (transfer_kB_eV * getParam<Real>("temperature"))) *
                 std::sqrt((1.0 - 0.5 * getParam<Real>("poisson_ratio")) /
                           (1.0 - getParam<Real>("poisson_ratio")) / libMesh::pi) /
                 (getParam<Real>("burgers_vector") * getParam<Real>("ca_crho"))),
    _minimum_size(getParam<unsigned int>("minimum_size")),
    _rho(getPostprocessorValue("dislocation_density_postprocessor"))
{
  if (getParam<Real>("D0") <= 0.0 || getParam<Real>("temperature") <= 0.0 ||
      getParam<Real>("burgers_vector") <= 0.0)
    mooseError(type(), " requires positive D0, temperature and burgers_vector.");
  if (getParam<Real>("bias") <= 0.0 || getParam<Real>("ca_crho") <= 0.0)
    mooseError(type(), " requires positive bias and ca_crho.");
  if (getParam<Real>("poisson_ratio") < 0.0 || getParam<Real>("poisson_ratio") >= 1.0)
    paramError("poisson_ratio", "must be in [0, 1).");
  if (_count < _minimum_size)
    paramError("minimum_size",
               "is ",
               _minimum_size,
               " but the variable has only ",
               _count,
               " components.");
}

Real
LoopNetworkTransferNodalKernel::rate(Real c) const
{
  return _coefficient * std::max(c, 0.0) * std::sqrt(std::max(_rho, 0.0));
}

void
LoopNetworkTransferNodalKernel::computeQpResidual(RealEigenVector & residual)
{
  residual.resize(_count);
  residual.setZero();
  const Real k = rate(_u[_qp](0));
  for (unsigned int i = _minimum_size - 1; i < _count; ++i)
    residual(i) = k * _u[_qp](i);
}

void
LoopNetworkTransferNodalKernel::computeQpJacobian()
{
  const Real k = rate(_u[_qp](0));
  const Real dk = _u[_qp](0) > 0.0 ? _coefficient * std::sqrt(std::max(_rho, 0.0)) : 0.0;
  for (unsigned int i = _minimum_size - 1; i < _count; ++i)
  {
    setJacobian(i, i, k);
    if (dk != 0.0)
      setJacobian(i, 0, dk * _u[_qp](i));
  }
}
