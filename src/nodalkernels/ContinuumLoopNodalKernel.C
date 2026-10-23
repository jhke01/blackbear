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

#include "FEProblemBase.h"
#include "Function.h"

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
      "to switch_size and a Fokker-Planck finite-volume continuum above it (Pokor et al. 2004, "
      "Section 4.1), optionally with the shrinking by and the loss of monomers of the opposite "
      "point-defect type.");

  params.addRequiredParam<Real>(
      "generation",
      "Total point-defect production rate [site fraction/s]. The fraction 1 - f2 - f3 - f4 of it "
      "is produced as monomers.");
  params.addParam<Real>("cascade_f2", 0.0, "Fraction of the production that is made as dimers.");
  params.addParam<Real>("cascade_f3", 0.0, "Fraction of the production that is made as trimers.");
  params.addParam<Real>(
      "cascade_f4", 0.0, "Fraction of the production that is made as clusters of 4.");
  params.addParam<Real>("sink", 0.0, "Linear sink coefficient of the monomer, k_s [1/s].");
  params.addParam<FunctionName>("sink_function",
                                "1",
                                "Function of time that multiplies sink, for a sink strength that "
                                "changes in time such as recovering dislocations.");
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

  params.addCoupledVar("opposite_coupled_var",
                       0.0,
                       "The other point-defect family's array variable, which must use the same "
                       "num_cluster_sizes, switch_size and cell widths. Optional: when absent the "
                       "opposite-species terms are inactive.");
  params.addParam<Real>("opposite_D0",
                        0.0,
                        "Diffusion prefactor of the opposite monomer [m^2/s] (needed with "
                        "opposite_coupled_var).");
  params.addParam<Real>("opposite_Em_eV",
                        0.0,
                        "Migration energy of the opposite monomer [eV] (needed with "
                        "opposite_coupled_var).");
  params.addParam<Real>("opposite_Z",
                        1.0,
                        "Large-loop limit of the bias of this family's loops for capturing the "
                        "opposite monomer.");
  params.addParam<Real>("opposite_Z_loop",
                        1.0,
                        "Small-loop coefficient of the bias of this family's loops for capturing "
                        "the opposite monomer.");
  params.addParam<Real>("opposite_alpha",
                        1.0,
                        "Exponent of the approach of that bias to opposite_Z with size.");
  params.addParam<Real>("opposite_burgers_vector",
                        0.0,
                        "Burgers vector of the opposite family's loops [m], which sets their "
                        "radius when they absorb this family's monomer. Defaults to burgers_vector "
                        "when zero.");
  return params;
}

ContinuumLoopNodalKernel::ContinuumLoopNodalKernel(const InputParameters & parameters)
  : ArrayNodalKernel(parameters),
    _generation(getParam<Real>("generation")),
    _cascade_fraction{1.0 - getParam<Real>("cascade_f2") - getParam<Real>("cascade_f3") -
                          getParam<Real>("cascade_f4"),
                      getParam<Real>("cascade_f2"),
                      getParam<Real>("cascade_f3"),
                      getParam<Real>("cascade_f4")},
    _base_sink(getParam<Real>("sink")),
    _sink_function(getFunction("sink_function")),
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
    _opposite_active(isCoupled("opposite_coupled_var")),
    _opposite_diffusivity(
        getParam<Real>("opposite_D0") *
        std::exp(-getParam<Real>("opposite_Em_eV") * continuum_eV / (continuum_kB * _temperature))),
    _opposite_burgers_vector(getParam<Real>("opposite_burgers_vector") != 0.0
                                 ? getParam<Real>("opposite_burgers_vector")
                                 : _burgers_vector),
    _opposite_z(getParam<Real>("opposite_Z")),
    _opposite_z_loop(getParam<Real>("opposite_Z_loop")),
    _opposite_alpha(getParam<Real>("opposite_alpha")),
    _opposite_array(coupledArrayValue("opposite_coupled_var")),
    _opposite_var_num(coupled("opposite_coupled_var")),
    _layout(getParam<unsigned int>("num_cluster_sizes"),
            getParam<unsigned int>("switch_size"),
            getParam<Real>("continuum_initial_width"),
            getParam<Real>("continuum_width_ratio")),
    _flux(_layout, getParam<bool>("use_limiter"))
{
  if (_temperature <= 0.0)
    paramError("temperature", "must be positive.");
  if (_cascade_fraction[0] < 0.0 || _cascade_fraction[1] < 0.0 || _cascade_fraction[2] < 0.0 ||
      _cascade_fraction[3] < 0.0)
    mooseError(type(), " requires cascade_f2, cascade_f3 and cascade_f4 >= 0 with a sum of at most 1.");
  if (_layout.switchSize() < 4)
    mooseError(type(), " requires switch_size >= 4 so that the cascade source is discrete.");
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

  if (_opposite_active)
  {
    if (_opposite_var_num == _var.number())
      paramError("opposite_coupled_var", "must be the other point-defect family's variable.");
    if (getParam<Real>("opposite_D0") <= 0.0)
      paramError("opposite_D0", "must be positive when opposite_coupled_var is given.");
    const unsigned int opposite_count = getArrayVar("opposite_coupled_var", 0)->count();
    if (opposite_count != _count)
      paramError("opposite_coupled_var",
                 "has ",
                 opposite_count,
                 " components but this variable has ",
                 _count,
                 "; the two families must use the same layout.");
    // The Jacobian has entries in every row of the monomer and in the shrink fluxes, outside the
    // diagonal blocks of the default pattern
    if (!_fe_problem.useHashTableMatrixAssembly())
      mooseError(type(), " with opposite_coupled_var requires Problem/use_hash_table_matrix_assembly = true");
  }
}

Real
ContinuumLoopNodalKernel::biasFactor(
    unsigned int n, Real burgers, Real z, Real z_loop, Real alpha_exponent) const
{
  Real bias = std::sqrt(burgers / (8.0 * libMesh::pi * _lattice_parameter)) * z_loop - z;
  bias *= std::pow(1.0 / static_cast<Real>(n), 0.5 * alpha_exponent);
  return bias + z;
}

Real
ContinuumLoopNodalKernel::currentSink() const
{
  return _base_sink * _sink_function.value(_t, *_current_node);
}

Real
ContinuumLoopNodalKernel::oppositeBeta(unsigned int n) const
{
  if (n == 0)
    return 0.0;
  const Real radius =
      std::sqrt(static_cast<Real>(n) * _atomic_volume / (libMesh::pi * _burgers_vector));
  return 2.0 * libMesh::pi * radius *
         biasFactor(n, _burgers_vector, _opposite_z, _opposite_z_loop, _opposite_alpha) *
         _opposite_diffusivity / _atomic_volume;
}

Real
ContinuumLoopNodalKernel::partnerBeta(unsigned int n) const
{
  if (n == 0)
    return 0.0;
  const Real radius =
      std::sqrt(static_cast<Real>(n) * _atomic_volume / (libMesh::pi * _opposite_burgers_vector));
  return 2.0 * libMesh::pi * radius * biasFactor(n, _burgers_vector, _z, _z_loop, _alpha_loop) *
         _diffusivity / _atomic_volume;
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
  if (_opposite_active)
  {
    _own_absorption.build(_layout, [this](unsigned int n) { return oppositeBeta(n); });
    _partner_absorption.build(_layout, [this](unsigned int n) { return partnerBeta(n); });
  }
  _have_coefficients = true;
}

void
ContinuumLoopNodalKernel::computeQpResidual(RealEigenVector & residual)
{
  ensureCoefficients();
  ClusterContinuumOpposite opposite{&_own_absorption, &_partner_absorption, nullptr};
  if (_opposite_active)
    opposite.partner_state = &_opposite_array[_qp](0);
  _flux.rate(&_u[_qp](0),
             _coefficients,
             _generation * _cascade_fraction[0],
             currentSink(),
             _rate,
             _opposite_active ? &opposite : nullptr);
  // A cluster of size n made in a cascade carries n point defects, so the production fraction f_n
  // makes f_n / n clusters
  for (const auto n : make_range(2u, 5u))
    _rate[n - 1] += _generation * _cascade_fraction[n - 1] / static_cast<Real>(n);
  residual.resize(_count);
  for (const auto i : make_range(_count))
    residual(i) = -_rate[i];
}

void
ContinuumLoopNodalKernel::computeQpJacobian()
{
  ensureCoefficients();
  ClusterContinuumOpposite opposite{&_own_absorption, &_partner_absorption, nullptr};
  if (_opposite_active)
    opposite.partner_state = &_opposite_array[_qp](0);
  _flux.jacobian(&_u[_qp](0),
                 _coefficients,
                 currentSink(),
                 _triplets,
                 _opposite_active ? &opposite : nullptr);
  for (const auto & t : _triplets)
    if (t.column < _count)
      setJacobian(t.row, t.column, t.value);
}

void
ContinuumLoopNodalKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (!_opposite_active || jvar != _opposite_var_num)
    return;

  ensureCoefficients();
  ClusterContinuumOpposite opposite{
      &_own_absorption, &_partner_absorption, &_opposite_array[_qp](0)};
  _flux.jacobian(&_u[_qp](0), _coefficients, currentSink(), _triplets, &opposite);
  for (const auto & t : _triplets)
    if (t.column >= _count)
      setJacobian(t.row, t.column - _count, t.value);
}
