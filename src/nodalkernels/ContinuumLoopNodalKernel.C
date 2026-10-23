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

#include "ContinuumLoopNodalKernel.h"

#include "libmesh/utility.h"

registerMooseObject("BlackBearApp", ContinuumLoopNodalKernel);

namespace
{
constexpr Real continuum_kB = 1.380649e-23;
constexpr Real continuum_eV = 1.602176634e-19;
} // namespace

InputParameters
ContinuumLoopNodalKernel::validParams()
{
  InputParameters params = ArrayNodalKernel::validParams();
  params.addClassDescription(
      "Absorption and emission of two-dimensional loops whose size distribution is discrete up "
      "to switch_size and a Fokker-Planck finite-volume continuum above it, for one family "
      "without opposite-species terms (Pokor et al. 2004, Section 4.1).");

  params.addRequiredParam<Real>("generation", "Monomer generation rate G_1 [1/s].");
  params.addParam<Real>("sink", 0.0, "Linear sink coefficient of the monomer, k_s [1/s].");
  params.addRequiredParam<Real>("temperature", "Temperature T [K].");
  params.addRequiredParam<Real>("atomic_volume", "Atomic volume V_at [m^3].");
  params.addRequiredParam<Real>("lattice_parameter", "Lattice parameter a0 [m].");
  params.addParam<Real>("burgers_vector",
                        0.0,
                        "Burgers vector b [m]. Defaults to lattice_parameter/sqrt(3) when zero.");
  params.addRequiredParam<Real>("D0", "Monomer diffusion prefactor [m^2/s].");
  params.addRequiredParam<Real>("Em_eV", "Monomer migration energy [eV].");
  params.addParam<Real>("Z", 1.0, "Large-loop limit of the self-absorption bias.");
  params.addParam<Real>("Z_loop", 1.0, "Small-loop coefficient of the self-absorption bias.");
  params.addParam<Real>("alpha", 1.0, "Exponent of the approach of the bias to Z with size.");
  params.addRequiredParam<Real>("formation_energy_eV", "Monomer formation energy E_f [eV].");
  params.addRequiredParam<Real>("energy_exponent", "Exponent p of the capillary law A n^p + C.");
  params.addRequiredParam<Real>("energy_coefficient_eV",
                                "Coefficient A of the capillary law [eV].");
  params.addRequiredParam<Real>("energy_constant_eV", "Constant C of the capillary law [eV].");

  params.addRequiredParam<unsigned int>("num_cluster_sizes",
                                        "Largest cluster size represented, N.");
  params.addRequiredParam<unsigned int>(
      "switch_size", "Largest cluster size that is an individual unknown; larger sizes are cells.");
  params.addParam<Real>("continuum_initial_width", 1.0, "Width of the first continuum cell.");
  params.addParam<Real>(
      "continuum_width_ratio", 1.04, "Growth factor of the cell width from one cell to the next.");
  params.addParam<bool>("use_limiter",
                        true,
                        "Reconstruct the cell values at a face with a minmod-limited slope "
                        "(second order, as in the SUNDIALS Pokor benchmark). If false the cell "
                        "averages are used (first-order upwind).");
  return params;
}

ContinuumLoopNodalKernel::ContinuumLoopNodalKernel(const InputParameters & parameters)
  : ArrayNodalKernel(parameters),
    _generation(getParam<Real>("generation")),
    _sink(getParam<Real>("sink")),
    _temperature(getParam<Real>("temperature")),
    _atomic_volume(getParam<Real>("atomic_volume")),
    _lattice_parameter(getParam<Real>("lattice_parameter")),
    _burgers_vector(getParam<Real>("burgers_vector") != 0.0
                        ? getParam<Real>("burgers_vector")
                        : getParam<Real>("lattice_parameter") / std::sqrt(3.0)),
    _diffusivity(getParam<Real>("D0") *
                 std::exp(-getParam<Real>("Em_eV") * continuum_eV / (continuum_kB * _temperature))),
    _z(getParam<Real>("Z")),
    _z_loop(getParam<Real>("Z_loop")),
    _alpha_loop(getParam<Real>("alpha")),
    _formation_energy_eV(getParam<Real>("formation_energy_eV")),
    _energy_exponent(getParam<Real>("energy_exponent")),
    _energy_coefficient_eV(getParam<Real>("energy_coefficient_eV")),
    _energy_constant_eV(getParam<Real>("energy_constant_eV")),
    _layout(getParam<unsigned int>("num_cluster_sizes"),
            getParam<unsigned int>("switch_size"),
            getParam<Real>("continuum_initial_width"),
            getParam<Real>("continuum_width_ratio")),
    _flux(_layout, getParam<bool>("use_limiter"))
{
  if (_temperature <= 0.0)
    paramError("temperature", "must be positive.");
  if (_atomic_volume <= 0.0 || _lattice_parameter <= 0.0 || _burgers_vector <= 0.0)
    mooseError(type(), " requires atomic_volume, lattice_parameter and burgers_vector > 0.");
  if (_energy_exponent <= 0.0 || _energy_coefficient_eV <= 0.0)
    mooseError(type(), " requires energy_exponent and energy_coefficient_eV > 0.");
  if (_count != _layout.componentCount())
    mooseError(type(),
               " expects ",
               _layout.componentCount(),
               " array components from the layout, but the variable has ",
               _count,
               ".");
}

Real
ContinuumLoopNodalKernel::beta(unsigned int n) const
{
  if (n == 0)
    return 0.0;
  const Real radius =
      std::sqrt(static_cast<Real>(n) * _atomic_volume / (libMesh::pi * _burgers_vector));
  // Z(n) = Z + [sqrt(b / (8 pi a0)) Z_loop - Z] n^(-alpha / 2): sqrt(b / (8 pi a0)) Z_loop at n =
  // 1, approaching Z for large loops
  Real bias = std::sqrt(_burgers_vector / (8.0 * libMesh::pi * _lattice_parameter)) * _z_loop - _z;
  bias *= std::pow(1.0 / static_cast<Real>(n), 0.5 * _alpha_loop);
  bias += _z;
  return 2.0 * libMesh::pi * radius * bias * _diffusivity / _atomic_volume;
}

Real
ContinuumLoopNodalKernel::alpha(unsigned int n) const
{
  if (n < 2)
    return 0.0;
  auto formation = [this](unsigned int m)
  {
    return _energy_coefficient_eV * std::pow(static_cast<Real>(m), _energy_exponent) +
           _energy_constant_eV;
  };
  const Real binding_eV = _formation_energy_eV - (formation(n) - formation(n - 1));
  // The exponent is clamped as in the SUNDIALS Pokor benchmark so that a very large binding energy
  // gives a zero emission rate instead of an underflow
  const Real exponent = -binding_eV * continuum_eV / (continuum_kB * _temperature);
  const Real boltzmann =
      exponent < -100.0 ? 0.0 : (exponent > 100.0 ? std::exp(100.0) : std::exp(exponent));
  return beta(n - 1) * boltzmann;
}

void
ContinuumLoopNodalKernel::ensureCoefficients()
{
  if (_have_coefficients)
    return;
  _coefficients.build(
      _layout,
      [this](unsigned int n) { return beta(n); },
      [this](unsigned int n) { return alpha(n); });
  _have_coefficients = true;
}

void
ContinuumLoopNodalKernel::computeQpResidual(RealEigenVector & residual)
{
  ensureCoefficients();
  _flux.rate(&_u[_qp](0), _coefficients, _generation, _sink, _rate);
  residual.resize(_count);
  for (const auto i : make_range(_count))
    residual(i) = -_rate[i];
}

void
ContinuumLoopNodalKernel::computeQpJacobian()
{
  ensureCoefficients();
  _flux.jacobian(&_u[_qp](0), _coefficients, _sink, _triplets);
  for (const auto & t : _triplets)
    setJacobian(t.row, t.column, t.value);
}
