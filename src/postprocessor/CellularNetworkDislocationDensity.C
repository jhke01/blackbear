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

#include "CellularNetworkDislocationDensity.h"

#include <cmath>

registerMooseObject("BlackBearApp", CellularNetworkDislocationDensity);

namespace
{
constexpr Real network_kB_eV = 8.617333262e-5;
}

InputParameters
CellularNetworkDislocationDensity::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addClassDescription(
      "Cellular network dislocation density of Rest (2006), Eqs 14-18, advanced at the end of "
      "each time step with the interstitial loop density at its steady state. Read by "
      "DislocationRecoveryNodalKernel through dislocation_density_postprocessor.");
  params.addRequiredParam<PostprocessorName>(
      "interstitial_concentration",
      "Postprocessor with the interstitial monomer concentration c_i [site fraction].");
  params.addRequiredParam<Real>("D0", "Interstitial diffusion prefactor D0 [m^2/s].");
  params.addRequiredParam<Real>("Em_eV", "Interstitial migration energy [eV].");
  params.addRequiredParam<Real>("temperature", "Temperature [K].");
  params.addRequiredParam<Real>("burgers_vector", "Burgers vector b [m] of the loops.");
  params.addRequiredParam<Real>("atomic_volume", "Atomic volume Omega [m^3].");
  params.addRequiredParam<Real>("rho0", "Initial network density [m^-2].");
  params.addRequiredParam<Real>(
      "bias",
      "Relative interstitial/vacancy bias Z_iv of Rest Eq 16. Material dependent; Rest's value for "
      "UO2 is 5e-4. Not known for steel.");
  params.addParam<Real>("ca_crho",
                        2.0,
                        "Product C_A C_rho of Rest Eq 10 and 14 (C_A = 3 for cubic cells, C_rho of "
                        "order one; 2 is the value of Alshannaq and Khafizov).");
  params.addParam<Real>("poisson_ratio", 0.3, "Poisson ratio nu in f(nu) = (1 - nu/2)/(1 - nu).");
  params.addParam<Real>("grain_size",
                        0.0,
                        "Grain size d_g [m] of the grain-boundary loss term 4 v_l rho/d_g. "
                        "0 turns the term off.");
  params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
  return params;
}

CellularNetworkDislocationDensity::CellularNetworkDislocationDensity(
    const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _c_i(getPostprocessorValue("interstitial_concentration")),
    _D0(getParam<Real>("D0")),
    _Em_eV(getParam<Real>("Em_eV")),
    _temperature(getParam<Real>("temperature")),
    _burgers_vector(getParam<Real>("burgers_vector")),
    _atomic_volume(getParam<Real>("atomic_volume")),
    _rho0(getParam<Real>("rho0")),
    _z_iv(getParam<Real>("bias")),
    _ca_crho(getParam<Real>("ca_crho")),
    _poisson_ratio(getParam<Real>("poisson_ratio")),
    _grain_size(getParam<Real>("grain_size")),
    _diffusivity(_D0 * std::exp(-_Em_eV / (network_kB_eV * _temperature))),
    _rho(declareRestartableData<Real>("network_density", _rho0)),
    _c_i_old(declareRestartableData<Real>("interstitial_concentration_old", 0.0))
{
  if (_D0 <= 0.0 || _temperature <= 0.0 || _burgers_vector <= 0.0 || _atomic_volume <= 0.0)
    mooseError(type(), " requires positive D0, temperature, burgers_vector and atomic_volume.");
  if (_rho0 <= 0.0)
    paramError("rho0", "must be positive.");
  if (_z_iv <= 0.0)
    paramError("bias", "must be positive.");
  if (_ca_crho <= 0.0)
    paramError("ca_crho", "must be positive.");
  if (_poisson_ratio < 0.0 || _poisson_ratio >= 1.0)
    paramError("poisson_ratio", "must be in [0, 1).");
  if (_grain_size < 0.0)
    paramError("grain_size", "must not be negative.");
}

Real
CellularNetworkDislocationDensity::advance(Real rho0,
                                           Real dt,
                                           Real c_i,
                                           Real diffusivity,
                                           Real burgers_vector,
                                           Real atomic_volume,
                                           Real z_iv,
                                           Real ca_crho,
                                           Real poisson_ratio,
                                           Real grain_size)
{
  if (c_i <= 0.0 || dt <= 0.0)
    return rho0;

  // Rest Eq 14 with Eq 18 and u = sqrt(rho):
  //   2 u du/dt = pi v k u - 4 v u^2/d_g u ... -> du/dt = (pi v k)/2 - (2 v/d_g) u - (C v/2) u^2,
  // with v = 2 Z_iv D_i c_i / b, k = [sqrt(2) b c_i/(4 pi Z_iv Omega^(5/3))]^(1/2),
  // C = (f/pi)^(1/2)/(C_A C_rho)
  const Real pi = libMesh::pi;
  const Real f = (1.0 - 0.5 * poisson_ratio) / (1.0 - poisson_ratio);
  const Real v = 2.0 * z_iv * diffusivity * c_i / burgers_vector;
  const Real k = std::sqrt(std::sqrt(2.0) * burgers_vector * c_i /
                           (4.0 * pi * z_iv * std::pow(atomic_volume, 5.0 / 3.0)));
  const Real C = std::sqrt(f / pi) / ca_crho;
  const Real a = 0.5 * C * v;                                     // coefficient of u^2
  const Real b = grain_size > 0.0 ? 2.0 * v / grain_size : 0.0;   // coefficient of u
  const Real c = 0.5 * pi * v * k;                                // source
  // du/dt = c - b u - a u^2 = -a (u - u_plus)(u - u_minus)
  const Real disc = std::sqrt(b * b + 4.0 * a * c);
  const Real u_plus = (-b + disc) / (2.0 * a);
  const Real u_minus = (-b - disc) / (2.0 * a);
  const Real u0 = std::sqrt(rho0);
  // w = (u - u_plus)/(u - u_minus) obeys dw/dt = -a (u_plus - u_minus) w
  const Real w0 = (u0 - u_plus) / (u0 - u_minus);
  const Real w = w0 * std::exp(-a * (u_plus - u_minus) * dt);
  const Real u = (u_plus - u_minus * w) / (1.0 - w);
  return u * u;
}

void
CellularNetworkDislocationDensity::execute()
{
  if (_t_step == 0)
  {
    _rho = _rho0;
    _c_i_old = _c_i;
    return;
  }
  const Real c_mean = 0.5 * (_c_i_old + _c_i);
  _rho = advance(_rho, _dt, c_mean, _diffusivity, _burgers_vector, _atomic_volume, _z_iv,
                 _ca_crho, _poisson_ratio, _grain_size);
  _c_i_old = _c_i;
}
