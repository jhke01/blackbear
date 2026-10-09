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

#include "GenericVacancyVoidGroupedNodalKernel.h"
#include "ClusterDynamicsSolverCheck.h"
#include "FEProblemBase.h"
#include "NonlinearSystemBase.h"

#include "libmesh/libmesh_common.h"

#include <cstddef>

namespace
{
constexpr Real void_kB = 1.380649e-23;
constexpr Real void_eV = 1.602176634e-19;
} // namespace

registerMooseObject("BlackBearApp", VacancyVoidGroupedNodalKernel);
registerMooseObject("BlackBearApp", ADVacancyVoidGroupedNodalKernel);

template <bool is_ad>
InputParameters
GenericVacancyVoidGroupedNodalKernelTempl<is_ad>::validParams()
{
  InputParameters params = GenericArrayNodalKernel<is_ad>::validParams();
  MooseEnum group_spacing("linear geometric", "linear");
  params.addClassDescription(
      "Grouped L0/L1 vacancy-void self-reaction nodal kernel: three-dimensional "
      "void growth/shrinkage by vacancy-monomer absorption and capillary-law "
      "thermal emission, including a finite-sink-strength correction to the "
      "absorption rate. Point-defect recombination with an interstitial "
      "variable is added separately by PointDefectRecombinationNodalKernel.");
  params.addParam<Real>("generation",
                        0.0,
                        "Vacancy monomer generation rate G_v1, e.g. dpa_rate * "
                        "cascade_efficiency for a monomer-only cascade spectrum.");
  params.addParam<Real>(
      "sink", 0.0, "Linear dislocation sink coefficient for vacancy monomer loss.");
  params.addRequiredParam<unsigned int>(
      "num_cluster_sizes",
      "Maximum physical void size (in vacancies) represented by the grouped distribution.");
  params.addParam<bool>(
      "use_grouping", true, "Whether to activate grouped-tail storage above explicit_max.");
  params.addParam<unsigned int>(
      "explicit_max", 100, "Maximum void size retained explicitly before grouped bins.");
  params.addParam<MooseEnum>(
      "group_spacing", group_spacing, "Grouped-tail spacing mode: linear or geometric.");
  params.addParam<unsigned int>(
      "group_linear_width", 25, "Fixed grouped-bin width when group_spacing = linear.");
  params.addParam<Real>("group_geometric_ratio",
                        1.35,
                        "Geometric grouped-bin growth factor when group_spacing = geometric.");
  params.addParam<bool>(
      "enforce_group_nonnegative",
      true,
      "Whether to limit grouped L0/L1 reconstruction to keep bin endpoint concentrations "
      "nonnegative (up to the tolerance factor).");
  params.addParam<Real>(
      "group_nonnegative_tolerance_factor",
      0.0,
      "Relaxation factor for grouped nonnegativity enforcement. Endpoint concentrations are "
      "allowed down to -factor * L0 before limiting.");
  params.addParam<Real>("temperature", 0.0, "Temperature T [K].");
  params.addParam<Real>("D0_v", 0.0, "Vacancy diffusion prefactor D0_v [m^2/s].");
  params.addParam<Real>("Em_v_eV", 0.0, "Vacancy migration energy Em_v [eV].");
  params.addParam<Real>("atomic_volume", 0.0, "Atomic volume V_at [m^3].");
  params.addParam<Real>(
      "lattice_parameter", 0.0, "Lattice parameter a0 [m] for the finite-sink-strength correction.");
  params.addParam<Real>(
      "formation_energy_eV", 0.0, "Vacancy monomer formation energy E_form_v [eV].");
  params.addParam<Real>(
      "void_energy_exponent", 2.0 / 3.0, "Exponent p_v in the capillary void formation energy law.");
  params.addParam<Real>("void_energy_coefficient_eV",
                        0.0,
                        "Coefficient A_v [eV] in the capillary void formation energy law "
                        "E_f(n) = A_v*n^p_v + C_v.");
  params.addParam<Real>("void_energy_constant_eV",
                        0.0,
                        "Constant C_v [eV] in the capillary void formation energy law "
                        "E_f(n) = A_v*n^p_v + C_v.");
  return params;
}

template <bool is_ad>
GenericVacancyVoidGroupedNodalKernelTempl<is_ad>::GenericVacancyVoidGroupedNodalKernelTempl(
    const InputParameters & parameters)
  : GenericArrayNodalKernel<is_ad>(parameters),
    _generation(this->template getParam<Real>("generation")),
    _sink(this->template getParam<Real>("sink")),
    _temperature(this->template getParam<Real>("temperature")),
    _D0(this->template getParam<Real>("D0_v")),
    _Q(this->template getParam<Real>("Em_v_eV") * void_eV),
    _atomic_volume(this->template getParam<Real>("atomic_volume")),
    _lattice_parameter(this->template getParam<Real>("lattice_parameter")),
    _formation_energy_eV(this->template getParam<Real>("formation_energy_eV")),
    _energy_exponent(this->template getParam<Real>("void_energy_exponent")),
    _energy_coefficient_eV(this->template getParam<Real>("void_energy_coefficient_eV")),
    _energy_constant_eV(this->template getParam<Real>("void_energy_constant_eV")),
    _num_cluster_sizes(this->template getParam<unsigned int>("num_cluster_sizes")),
    _use_grouping(this->template getParam<bool>("use_grouping")),
    _explicit_max(this->template getParam<unsigned int>("explicit_max")),
    _group_spacing(this->template getParam<MooseEnum>("group_spacing")
                       .template getEnum<ClusterGroupingSpacing>()),
    _group_linear_width(this->template getParam<unsigned int>("group_linear_width")),
    _group_geometric_ratio(this->template getParam<Real>("group_geometric_ratio")),
    _enforce_group_nonnegative(this->template getParam<bool>("enforce_group_nonnegative")),
    _group_nonnegative_tolerance_factor(
        this->template getParam<Real>("group_nonnegative_tolerance_factor")),
    _layout(_num_cluster_sizes,
            _explicit_max,
            _use_grouping,
            _group_spacing,
            _group_linear_width,
            _group_geometric_ratio)
{
  // Grouped and explicit components are coupled through the monomer and the inter-bin fluxes,
  // so the Jacobian contains entries outside the default diagonal-block sparsity pattern.
  if (!this->_fe_problem.useHashTableMatrixAssembly())
    mooseError(
        "VacancyVoidGroupedNodalKernel requires Problem/use_hash_table_matrix_assembly = true");
  if (this->_tid == 0)
    ClusterDynamics::checkSolverSetup(*this,
                                      this->_fe_problem,
                                      this->_var.sys().number(),
                                      this->_var.number(),
                                      this->_var.count());

  if (_group_nonnegative_tolerance_factor < 0.0)
    mooseError("VacancyVoidGroupedNodalKernel requires group_nonnegative_tolerance_factor >= 0.");
  if (_temperature <= 0.0)
    mooseError("VacancyVoidGroupedNodalKernel requires temperature > 0.");
  if (_D0 <= 0.0)
    mooseError("VacancyVoidGroupedNodalKernel requires D0_v > 0.");
  if (_Q <= 0.0)
    mooseError("VacancyVoidGroupedNodalKernel requires Em_v_eV > 0.");
  if (_atomic_volume <= 0.0)
    mooseError("VacancyVoidGroupedNodalKernel requires atomic_volume > 0.");
  if (_lattice_parameter <= 0.0)
    mooseError("VacancyVoidGroupedNodalKernel requires lattice_parameter > 0.");
  if (_formation_energy_eV <= 0.0)
    mooseError("VacancyVoidGroupedNodalKernel requires formation_energy_eV > 0.");
  if (_energy_exponent <= 0.0)
    mooseError("VacancyVoidGroupedNodalKernel requires void_energy_exponent > 0.");
  if (_energy_coefficient_eV <= 0.0)
    mooseError("VacancyVoidGroupedNodalKernel requires void_energy_coefficient_eV > 0.");
}

template <bool is_ad>
Real
GenericVacancyVoidGroupedNodalKernelTempl<is_ad>::radius(unsigned int n) const
{
  return std::cbrt(3.0 * static_cast<Real>(n) * _atomic_volume / (4.0 * libMesh::pi));
}

template <bool is_ad>
Real
GenericVacancyVoidGroupedNodalKernelTempl<is_ad>::formationEnergy(unsigned int n) const
{
  return _energy_coefficient_eV * std::pow(static_cast<Real>(n), _energy_exponent) +
         _energy_constant_eV;
}

template <bool is_ad>
Real
GenericVacancyVoidGroupedNodalKernelTempl<is_ad>::monomerDiffusivity() const
{
  return _D0 * std::exp(-_Q / (void_kB * _temperature));
}

template <bool is_ad>
Real
GenericVacancyVoidGroupedNodalKernelTempl<is_ad>::beta(unsigned int n) const
{
  if (n < _beta_cache.size())
    return _beta_cache[n];

  const Real r = radius(n);
  return 4.0 * libMesh::pi * r * monomerDiffusivity() /
         ((1.0 + _lattice_parameter / r) * _atomic_volume);
}

template <bool is_ad>
Real
GenericVacancyVoidGroupedNodalKernelTempl<is_ad>::alpha(unsigned int n) const
{
  if (n < _alpha_cache.size())
    return _alpha_cache[n];

  if (n <= 1)
    mooseError("alpha(n) is only defined for void sizes n >= 2.");

  const Real binding_energy_eV =
      _formation_energy_eV - (formationEnergy(n) - formationEnergy(n - 1));
  return beta(n - 1) * std::exp(-binding_energy_eV * void_eV / (void_kB * _temperature));
}

template <bool is_ad>
void
GenericVacancyVoidGroupedNodalKernelTempl<is_ad>::ensureCoefficientCache() const
{
  if (_cache_size == _num_cluster_sizes)
    return;

  // Compute into local vectors first: beta(n)/alpha(n) short-circuit to a cache lookup once
  // _beta_cache/_alpha_cache have the right size, so filling them in place while they are already
  // sized would just read back the zeros they were assigned.
  std::vector<Real> beta_values(_num_cluster_sizes + 1, 0.0);
  std::vector<Real> alpha_values(_num_cluster_sizes + 1, 0.0);
  const Real diffusivity = monomerDiffusivity();

  for (unsigned int n = 1; n <= _num_cluster_sizes; ++n)
  {
    const Real r = radius(n);
    beta_values[n] =
        4.0 * libMesh::pi * r * diffusivity / ((1.0 + _lattice_parameter / r) * _atomic_volume);
  }
  for (unsigned int n = 2; n <= _num_cluster_sizes; ++n)
  {
    const Real binding_energy_eV =
        _formation_energy_eV - (formationEnergy(n) - formationEnergy(n - 1));
    alpha_values[n] =
        beta_values[n - 1] * std::exp(-binding_energy_eV * void_eV / (void_kB * _temperature));
  }

  _beta_cache = std::move(beta_values);
  _alpha_cache = std::move(alpha_values);
  _cache_size = _num_cluster_sizes;
}

template <bool is_ad>
template <typename Scalar>
void
GenericVacancyVoidGroupedNodalKernelTempl<is_ad>::fillConcentrationCache(
    const std::vector<Scalar> & state, std::vector<Scalar> & concentration) const
{
  concentration.assign(_num_cluster_sizes + 1, Scalar(0.0));
  for (unsigned int n = 1; n <= _num_cluster_sizes; ++n)
    concentration[n] = _layout.concentration(
        state, n, _enforce_group_nonnegative, _group_nonnegative_tolerance_factor);
}

template <bool is_ad>
template <typename Scalar>
void
GenericVacancyVoidGroupedNodalKernelTempl<is_ad>::fillFluxCache(
    const std::vector<Scalar> & concentration, std::vector<Scalar> & flux) const
{
  flux.assign(_num_cluster_sizes + 1, Scalar(0.0));
  const Scalar c1 = concentration[1];
  for (unsigned int n = 1; n < _num_cluster_sizes; ++n)
    flux[n] = beta(n) * c1 * concentration[n] - alpha(n + 1) * concentration[n + 1];
}

template <bool is_ad>
template <typename Scalar>
void
GenericVacancyVoidGroupedNodalKernelTempl<is_ad>::computeGroupedResidualFromState(
    const std::vector<Scalar> & state, std::vector<Scalar> & residual_values) const
{
  std::vector<Scalar> concentration;
  fillConcentrationCache(state, concentration);
  std::vector<Scalar> flux;
  fillFluxCache(concentration, flux);

  residual_values.assign(_layout.componentCount(), Scalar(0.0));

  Scalar monomer_flux = 2.0 * flux[1];
  for (unsigned int n = 2; n < _num_cluster_sizes; ++n)
    monomer_flux += flux[n];
  residual_values[_layout.monomerComponent()] =
      -(_generation - _sink * concentration[1] - monomer_flux);

  for (unsigned int n = 2; n <= _layout.explicitMax(); ++n)
    residual_values[_layout.explicitComponent(n)] =
        -((n < _num_cluster_sizes ? flux[n - 1] - flux[n] : flux[n - 1]));

  for (const auto & bin : _layout.groups())
  {
    const auto j_left = flux[bin.start - 1];
    const auto j_right = bin.end < _num_cluster_sizes ? flux[bin.end] : Scalar(0.0);

    residual_values[bin.l0_component] = -((j_left - j_right) / bin.width);

    if (bin.sigma2 > 0.0 && bin.width > 1.0)
    {
      Scalar j_avg = 0.0;
      for (unsigned int n = bin.start; n < bin.end; ++n)
        j_avg += flux[n];
      j_avg /= (bin.width - 1.0);

      residual_values[bin.l1_component] =
          (j_left - 2.0 * j_avg + j_right) * (bin.width - 1.0) / (2.0 * bin.sigma2 * bin.width);
    }
    else
      residual_values[bin.l1_component] = 0.0;
  }
}

template <bool is_ad>
void
GenericVacancyVoidGroupedNodalKernelTempl<is_ad>::jacobianSetup()
{
  GenericArrayNodalKernel<is_ad>::jacobianSetup();

  // The preconditioner is configured from the PETSc options only once the solve starts, so the
  // factorization ordering is checked before the first Jacobian is factored
  if (_checked_factor_ordering || this->_tid != 0)
    return;
  if (auto * nl = dynamic_cast<NonlinearSystemBase *>(&this->_sys))
    _checked_factor_ordering = ClusterDynamics::checkFactorOrdering(*this, *nl);
}

template <bool is_ad>
void
GenericVacancyVoidGroupedNodalKernelTempl<is_ad>::computeQpResidual(
    GenericRealEigenVector<is_ad> & residual)
{
  const auto n_comp = _u[_qp].size();
  if (n_comp != _layout.componentCount())
    mooseError("VacancyVoidGroupedNodalKernel expected ",
               _layout.componentCount(),
               " array components from the grouping layout, but the variable has ",
               n_comp,
               ".");

  ensureCoefficientCache();
  std::vector<GenericReal<is_ad>> state(n_comp, GenericReal<is_ad>(0.0));
  for (unsigned int i = 0; i < n_comp; ++i)
    state[i] = _u[_qp](i);

  std::vector<GenericReal<is_ad>> residual_values;
  computeGroupedResidualFromState(state, residual_values);
  residual.resize(n_comp);
  for (unsigned int i = 0; i < n_comp; ++i)
    residual(i) = residual_values[i];
}

// Non-AD Jacobian: every residual row is a linear combination of the fluxes
// J_n = beta(n)*C_1*C_n - alpha(n+1)*C_{n+1}, and every reconstructed concentration is linear in
// the state, C_n = a0(n)*L0 + a1(n)*L1 for grouped sizes (with the limiter coefficients frozen).
// The chain rule through the fluxes gives the exact intra-variable Jacobian in O(N) work.
template <>
void
GenericVacancyVoidGroupedNodalKernelTempl<false>::computeQpJacobian()
{
  const auto n_comp = static_cast<unsigned int>(_u[_qp].size());
  if (n_comp != _layout.componentCount())
    mooseError("VacancyVoidGroupedNodalKernel expected ",
               _layout.componentCount(),
               " array components from the grouping layout, but the variable has ",
               n_comp,
               ".");

  ensureCoefficientCache();

  const RealEigenVector & state = _u[_qp];
  std::vector<Real> concentration(_num_cluster_sizes + 1, 0.0);
  std::vector<Real> a0(_num_cluster_sizes + 1, 1.0);
  std::vector<Real> a1(_num_cluster_sizes + 1, 0.0);
  for (unsigned int n = 1; n <= _num_cluster_sizes; ++n)
  {
    if (_layout.isExplicitSize(n))
    {
      concentration[n] = state(_layout.explicitComponent(n));
      continue;
    }
    const auto & bin = _layout.groupForSize(n);
    const Real l0 = state(bin.l0_component);
    const Real l1 = state(bin.l1_component);
    ClusterGroupingLayout::reconstructionCoefficients(bin,
                                                      n,
                                                      l0,
                                                      l1,
                                                      _enforce_group_nonnegative,
                                                      _group_nonnegative_tolerance_factor,
                                                      a0[n],
                                                      a1[n]);
    concentration[n] = a0[n] * l0 + a1[n] * l1;
  }
  const Real c1 = concentration[1];

  // Accumulate one Jacobian row at a time, visiting only the columns it touches
  std::vector<Real> row_values(n_comp, 0.0);
  std::vector<bool> row_touched(n_comp, false);
  std::vector<unsigned int> row_columns;
  auto addEntry = [&](unsigned int column, Real value)
  {
    if (!row_touched[column])
    {
      row_touched[column] = true;
      row_columns.push_back(column);
    }
    row_values[column] += value;
  };
  auto addConcentrationDerivative = [&](unsigned int n, Real value)
  {
    if (_layout.isExplicitSize(n))
      addEntry(_layout.explicitComponent(n), value);
    else
    {
      const auto & bin = _layout.groupForSize(n);
      addEntry(bin.l0_component, value * a0[n]);
      addEntry(bin.l1_component, value * a1[n]);
    }
  };
  auto addFluxDerivative = [&](unsigned int n, Real weight)
  {
    // J_0 and J_N vanish at the ends of the truncated size space
    if (n == 0 || n >= _num_cluster_sizes)
      return;
    addConcentrationDerivative(1, weight * beta(n) * concentration[n]);
    addConcentrationDerivative(n, weight * beta(n) * c1);
    addConcentrationDerivative(n + 1, -weight * alpha(n + 1));
  };
  auto setRow = [&](unsigned int row)
  {
    for (const auto column : row_columns)
    {
      setJacobian(row, column, row_values[column]);
      row_values[column] = 0.0;
      row_touched[column] = false;
    }
    row_columns.clear();
  };

  // Monomer: F = -(G_1 - k_s*C_1 - 2*J_1 - sum_{n=2}^{N-1} J_n)
  addEntry(_layout.monomerComponent(), _sink);
  for (unsigned int n = 1; n < _num_cluster_sizes; ++n)
    addFluxDerivative(n, n == 1 ? 2.0 : 1.0);
  setRow(_layout.monomerComponent());

  // Explicit sizes: F = -(J_{n-1} - J_n)
  for (unsigned int n = 2; n <= _layout.explicitMax(); ++n)
  {
    addFluxDerivative(n - 1, -1.0);
    addFluxDerivative(n, 1.0);
    setRow(_layout.explicitComponent(n));
  }

  for (const auto & bin : _layout.groups())
  {
    // L0: F = -(J_left - J_right) / width
    addFluxDerivative(bin.start - 1, -1.0 / bin.width);
    addFluxDerivative(bin.end, 1.0 / bin.width);
    setRow(bin.l0_component);

    // L1: F = (J_left - 2*J_avg + J_right) * (width - 1) / (2*sigma2*width)
    if (bin.sigma2 > 0.0 && bin.width > 1.0)
    {
      const Real factor = (bin.width - 1.0) / (2.0 * bin.sigma2 * bin.width);
      addFluxDerivative(bin.start - 1, factor);
      for (unsigned int n = bin.start; n < bin.end; ++n)
        addFluxDerivative(n, -2.0 * factor / (bin.width - 1.0));
      addFluxDerivative(bin.end, factor);
      setRow(bin.l1_component);
    }
  }
}

template <>
void
GenericVacancyVoidGroupedNodalKernelTempl<true>::computeQpJacobian()
{
  mooseError("Internal error: computeQpJacobian should never be called for the AD version");
}

template class GenericVacancyVoidGroupedNodalKernelTempl<false>;
template class GenericVacancyVoidGroupedNodalKernelTempl<true>;
