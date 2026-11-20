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

#include "ClusterDynamicsSolverCheck.h"
#include "FEProblemBase.h"

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
      "Absorption and emission of a family of two-dimensional loops or three-dimensional voids "
      "whose size distribution is discrete up to switch_size and a Fokker-Planck finite-volume "
      "continuum above it (Pokor et al. 2004, Section 4.1), optionally with the shrinking by and "
      "the loss of monomers of the opposite point-defect type.");

  params.addParam<MooseEnum>(
      "geometry",
      MooseEnum("loop void", "loop"),
      "Cluster geometry of this family: 'loop' (two-dimensional, radius ~ sqrt(n), dislocation "
      "bias) or 'void' (three-dimensional, radius ~ n^(1/3), finite-sink-strength absorption "
      "4 pi r D / ((1 + a0 / r) V_at)). A void ignores Z_loop and alpha, and Z is then only the "
      "large-loop limit of the opposite loops' bias for capturing this family's monomer.");
  params.addParam<MooseEnum>(
      "opposite_geometry",
      MooseEnum("loop void", "loop"),
      "Cluster geometry of the opposite family (opposite_coupled_var): sets the coefficient at "
      "which its clusters absorb this family's monomer.");
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
  params.addParam<Real>("partner_Z_loop",
                        "Small-loop coefficient of the bias of the opposite family's loops for "
                        "capturing THIS family's monomer. Defaults to Z_loop. The large-loop limit "
                        "is Z.");
  params.addParam<Real>("partner_alpha",
                        "Exponent of the approach of that bias to Z with size. Defaults to alpha.");
  params.addParam<Real>("opposite_burgers_vector",
                        0.0,
                        "Burgers vector of the opposite family's loops [m], which sets their "
                        "radius when they absorb this family's monomer. Defaults to burgers_vector "
                        "when zero.");
  params.addParam<bool>(
      "assemble_opposite_jacobian_in_own_block",
      false,
      "EXPERIMENTAL. If true, the cross-variable Jacobian block against opposite_coupled_var is "
      "written from this kernel's own-variable Jacobian call, using the other variable's dof "
      "indices directly, instead of through computeQpOffDiagJacobian. Use it with an SMP "
      "preconditioner that does NOT list the two variables in off_diag_row/off_diag_column: MOOSE "
      "then sees no coupling between different variables and does not allocate a local block for "
      "every pair of components of the whole system, while the assembled matrix still has the full "
      "cross-species coupling.");
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
    _is_void(getParam<MooseEnum>("geometry") == "void"),
    _opposite_is_void(getParam<MooseEnum>("opposite_geometry") == "void"),
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
    _partner_z_loop(isParamValid("partner_Z_loop") ? getParam<Real>("partner_Z_loop") : _z_loop),
    _partner_alpha(isParamValid("partner_alpha") ? getParam<Real>("partner_alpha") : _alpha_loop),
    _opposite_array(coupledArrayValue("opposite_coupled_var")),
    _opposite_var_num(coupled("opposite_coupled_var")),
    _assemble_opposite_in_own_block(getParam<bool>("assemble_opposite_jacobian_in_own_block")),
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

  // The monomer row and column are dense, so the Jacobian has entries outside the default pattern
  if (!_fe_problem.useHashTableMatrixAssembly())
    mooseError(type(), " requires Problem/use_hash_table_matrix_assembly = true");
  if (_tid == 0)
    ClusterDynamics::checkSolverSetup(
        *this, _fe_problem, _var.sys().number(), _var.number(), _var.count());

  if (_opposite_active)
  {
    if (_opposite_var_num == _var.number())
      paramError("opposite_coupled_var", "must be the other point-defect family's variable.");
    if (_tid == 0)
      ClusterDynamics::checkOppositeCoupling(*this,
                                             _fe_problem,
                                             _var.sys().number(),
                                             _var.number(),
                                             _opposite_var_num,
                                             _assemble_opposite_in_own_block);
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
ContinuumLoopNodalKernel::voidBeta(unsigned int n, Real diffusivity) const
{
  const Real radius = std::cbrt(3.0 * static_cast<Real>(n) * _atomic_volume / (4.0 * libMesh::pi));
  return 4.0 * libMesh::pi * radius * diffusivity /
         ((1.0 + _lattice_parameter / radius) * _atomic_volume);
}

Real
ContinuumLoopNodalKernel::oppositeBeta(unsigned int n) const
{
  if (n == 0)
    return 0.0;
  if (_is_void)
    return voidBeta(n, _opposite_diffusivity);
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
  if (_opposite_is_void)
    return voidBeta(n, _diffusivity);
  const Real radius =
      std::sqrt(static_cast<Real>(n) * _atomic_volume / (libMesh::pi * _opposite_burgers_vector));
  // The opposite family's own Burgers vector sets both the radius and the small-loop scale of the bias
  return 2.0 * libMesh::pi * radius *
         biasFactor(n, _opposite_burgers_vector, _z, _partner_z_loop, _partner_alpha) *
         _diffusivity / _atomic_volume;
}

Real
ContinuumLoopNodalKernel::beta(unsigned int n) const
{
  if (n == 0)
    return 0.0;
  if (_is_void)
    return voidBeta(n, _diffusivity);
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
             _generation,
             _sink,
             _rate,
             _opposite_active ? &opposite : nullptr);
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
                 _sink,
                 _triplets,
                 _opposite_active ? &opposite : nullptr);
  for (const auto & t : _triplets)
    if (t.column < _count)
      setJacobian(t.row, t.column, t.value);

  // Experimental: the cross-variable block, from the same triplets (the columns at or above
  // _count belong to the opposite family), written with the other variable's dof indices so that
  // no off_diag_row/off_diag_column pair is needed in the preconditioner
  if (_assemble_opposite_in_own_block && _opposite_active)
  {
    const auto & opposite_dofs = getVariable(_opposite_var_num).dofIndices();
    for (const auto & t : _triplets)
      if (t.column >= _count)
        addJacobianElement(_assembly,
                           t.value,
                           _var.dofIndices()[t.row],
                           opposite_dofs[t.column - _count],
                           _var.arrayScalingFactor()[t.row]);
  }
}

void
ContinuumLoopNodalKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (!_opposite_active || jvar != _opposite_var_num)
    return;
  // In the experimental mode the block is written from the own-variable call instead
  if (_assemble_opposite_in_own_block)
    return;

  ensureCoefficients();
  ClusterContinuumOpposite opposite{
      &_own_absorption, &_partner_absorption, &_opposite_array[_qp](0)};
  _flux.jacobian(&_u[_qp](0), _coefficients, _sink, _triplets, &opposite);
  for (const auto & t : _triplets)
    if (t.column >= _count)
      setJacobian(t.row, t.column - _count, t.value);
}
