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

#include "GenericInterstitialLoopGroupedNodalKernel.h"
#include "ClusterDynamicsSolverCheck.h"
#include "FEProblemBase.h"
#include "NonlinearSystemBase.h"

#include "libmesh/libmesh_common.h"

#include <cstddef>

namespace
{
constexpr Real loop_kB = 1.380649e-23;
constexpr Real loop_eV = 1.602176634e-19;
} // namespace

registerMooseObject("BlackBearApp", InterstitialLoopGroupedNodalKernel);
registerMooseObject("BlackBearApp", ADInterstitialLoopGroupedNodalKernel);

template <bool is_ad>
InputParameters
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::validParams()
{
  InputParameters params = GenericArrayNodalKernel<is_ad>::validParams();
  MooseEnum group_spacing("linear geometric", "linear");
  params.addClassDescription(
      "Grouped L0/L1 interstitial-loop self-reaction nodal kernel: two-dimensional "
      "loop growth/shrinkage by interstitial-monomer absorption and capillary-law "
      "thermal emission, including a size-dependent dislocation-loop bias factor "
      "in the absorption rate. Point-defect recombination with a vacancy variable "
      "is added separately by PointDefectRecombinationNodalKernel.");
  params.addParam<Real>("generation",
                        0.0,
                        "Interstitial monomer generation rate G_i1, e.g. dpa_rate * "
                        "cascade_efficiency for a monomer-only cascade spectrum.");
  params.addParam<Real>(
      "sink", 0.0, "Linear dislocation sink coefficient for interstitial monomer loss.");
  params.addRequiredParam<unsigned int>(
      "num_cluster_sizes",
      "Maximum physical loop size (in interstitials) represented by the grouped distribution.");
  params.addParam<bool>(
      "use_grouping", true, "Whether to activate grouped-tail storage above explicit_max.");
  params.addParam<unsigned int>(
      "explicit_max", 100, "Maximum loop size retained explicitly before grouped bins.");
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
  params.addParam<Real>("D0_i", 0.0, "Interstitial diffusion prefactor D0_i [m^2/s].");
  params.addParam<Real>("Em_i_eV", 0.0, "Interstitial migration energy Em_i [eV].");
  params.addParam<Real>("atomic_volume", 0.0, "Atomic volume V_at [m^3].");
  params.addParam<Real>(
      "lattice_parameter", 0.0, "Lattice parameter a0 [m] used to derive the Burgers vector.");
  params.addParam<Real>("burgers_vector",
                        0.0,
                        "Burgers vector b [m] for the loop radius and bias factor. Defaults to "
                        "lattice_parameter/sqrt(3) when left at zero.");
  params.addParam<Real>("z_i", 1.0, "Baseline dislocation bias Z_i for small loops.");
  params.addParam<Real>(
      "z_loop_i", 1.0, "Large-loop dislocation bias Z_loop_i approached as n -> infinity.");
  params.addParam<Real>(
      "alpha_loop_i", 1.0, "Exponent controlling how fast Z_i(n) approaches Z_loop_i with size.");
  params.addParam<Real>(
      "formation_energy_eV", 0.0, "Interstitial monomer formation energy E_form_i [eV].");
  params.addParam<Real>(
      "loop_energy_exponent", 2.0 / 3.0, "Exponent p_i in the capillary loop formation energy law.");
  params.addParam<Real>("loop_energy_coefficient_eV",
                        0.0,
                        "Coefficient A_i [eV] in the capillary loop formation energy law "
                        "E_f(n) = A_i*n^p_i + C_i.");
  params.addParam<Real>("loop_energy_constant_eV",
                        0.0,
                        "Constant C_i [eV] in the capillary loop formation energy law "
                        "E_f(n) = A_i*n^p_i + C_i.");

  params.addParam<Real>(
      "D0_v", 0.0, "Vacancy diffusion prefactor D0_v [m^2/s] (opposite species), for this "
      "loop's own absorption of a vacancy monomer.");
  params.addParam<Real>(
      "Em_v_eV", 0.0, "Vacancy migration energy Em_v [eV] (opposite species).");
  params.addParam<Real>(
      "z_v", 1.0, "Baseline dislocation bias Z_v for small loops absorbing vacancies.");
  params.addParam<Real>(
      "z_loop_v",
      1.0,
      "Large-loop dislocation bias Z_loop_v approached as n -> infinity, for loops absorbing "
      "vacancies.");
  params.addParam<Real>(
      "alpha_loop_v",
      1.0,
      "Exponent controlling how fast the loop vacancy-absorption bias approaches Z_loop_v with "
      "size.");

  params.addCoupledVar(
      "opposite_coupled_var",
      0.0,
      "The vacancy-void grouped array variable. Optional: when not provided, the "
      "opposite-species absorption terms (Piece A and Piece B) are inactive and this kernel "
      "reduces exactly to its self-reaction-only behavior.");
  params.addParam<unsigned int>(
      "opposite_num_cluster_sizes",
      1,
      "Maximum physical void size represented by the coupled grouped distribution. Only "
      "meaningful when opposite_coupled_var is provided.");
  params.addParam<bool>(
      "opposite_use_grouping", true, "Whether the coupled void array uses grouped-tail storage.");
  params.addParam<unsigned int>(
      "opposite_explicit_max", 100, "Maximum void size stored explicitly in the coupled array.");
  params.addParam<MooseEnum>(
      "opposite_group_spacing", group_spacing, "Grouped-tail spacing mode for the coupled array.");
  params.addParam<unsigned int>("opposite_group_linear_width",
                                25,
                                "Fixed grouped-bin width for the coupled array when "
                                "opposite_group_spacing = linear.");
  params.addParam<Real>("opposite_group_geometric_ratio",
                        1.35,
                        "Geometric grouped-bin growth factor for the coupled array when "
                        "opposite_group_spacing = geometric.");
  params.addParam<bool>("opposite_enforce_group_nonnegative",
                        true,
                        "Whether the coupled array's grouped reconstruction enforces "
                        "nonnegative bin endpoints.");
  params.addParam<Real>("opposite_group_nonnegative_tolerance_factor",
                        0.0,
                        "Relaxation factor for the coupled array's grouped nonnegativity "
                        "enforcement.");
  return params;
}

template <bool is_ad>
GenericInterstitialLoopGroupedNodalKernelTempl<
    is_ad>::GenericInterstitialLoopGroupedNodalKernelTempl(const InputParameters & parameters)
  : GenericArrayNodalKernel<is_ad>(parameters),
    _generation(this->template getParam<Real>("generation")),
    _sink(this->template getParam<Real>("sink")),
    _temperature(this->template getParam<Real>("temperature")),
    _D0(this->template getParam<Real>("D0_i")),
    _Q(this->template getParam<Real>("Em_i_eV") * loop_eV),
    _atomic_volume(this->template getParam<Real>("atomic_volume")),
    _lattice_parameter(this->template getParam<Real>("lattice_parameter")),
    _burgers_vector(this->template getParam<Real>("burgers_vector") != 0.0
                        ? this->template getParam<Real>("burgers_vector")
                        : this->template getParam<Real>("lattice_parameter") / std::sqrt(3.0)),
    _z_i(this->template getParam<Real>("z_i")),
    _z_loop_i(this->template getParam<Real>("z_loop_i")),
    _alpha_loop_i(this->template getParam<Real>("alpha_loop_i")),
    _formation_energy_eV(this->template getParam<Real>("formation_energy_eV")),
    _energy_exponent(this->template getParam<Real>("loop_energy_exponent")),
    _energy_coefficient_eV(this->template getParam<Real>("loop_energy_coefficient_eV")),
    _energy_constant_eV(this->template getParam<Real>("loop_energy_constant_eV")),
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
            _group_geometric_ratio),
    _opposite_D0(this->template getParam<Real>("D0_v")),
    _opposite_Q(this->template getParam<Real>("Em_v_eV") * loop_eV),
    _z_v(this->template getParam<Real>("z_v")),
    _z_loop_v(this->template getParam<Real>("z_loop_v")),
    _alpha_loop_v(this->template getParam<Real>("alpha_loop_v")),
    _opposite_num_cluster_sizes(this->template getParam<unsigned int>("opposite_num_cluster_sizes")),
    _opposite_use_grouping(this->template getParam<bool>("opposite_use_grouping")),
    _opposite_explicit_max(this->template getParam<unsigned int>("opposite_explicit_max")),
    _opposite_group_spacing(this->template getParam<MooseEnum>("opposite_group_spacing")
                                .template getEnum<ClusterGroupingSpacing>()),
    _opposite_group_linear_width(
        this->template getParam<unsigned int>("opposite_group_linear_width")),
    _opposite_group_geometric_ratio(
        this->template getParam<Real>("opposite_group_geometric_ratio")),
    _opposite_enforce_group_nonnegative(
        this->template getParam<bool>("opposite_enforce_group_nonnegative")),
    _opposite_group_nonnegative_tolerance_factor(
        this->template getParam<Real>("opposite_group_nonnegative_tolerance_factor")),
    _opposite_layout(_opposite_num_cluster_sizes,
                     _opposite_explicit_max,
                     _opposite_use_grouping,
                     _opposite_group_spacing,
                     _opposite_group_linear_width,
                     _opposite_group_geometric_ratio),
    _opposite_array(this->coupledArrayValue("opposite_coupled_var")),
    _opposite_var_num(this->coupled("opposite_coupled_var")),
    _opposite_coupling_active(this->isCoupled("opposite_coupled_var"))
{
  // Grouped and explicit components are coupled through the monomer and the inter-bin fluxes,
  // so the Jacobian contains entries outside the default diagonal-block sparsity pattern.
  if (!this->_fe_problem.useHashTableMatrixAssembly())
    mooseError("InterstitialLoopGroupedNodalKernel requires Problem/use_hash_table_matrix_"
               "assembly = true");
  if (this->_tid == 0)
    ClusterDynamics::checkSolverSetup(*this,
                                      this->_fe_problem,
                                      this->_var.sys().number(),
                                      this->_var.number(),
                                      this->_var.count());

  if (_group_nonnegative_tolerance_factor < 0.0)
    mooseError(
        "InterstitialLoopGroupedNodalKernel requires group_nonnegative_tolerance_factor >= 0.");
  if (_temperature <= 0.0)
    mooseError("InterstitialLoopGroupedNodalKernel requires temperature > 0.");
  if (_D0 <= 0.0)
    mooseError("InterstitialLoopGroupedNodalKernel requires D0_i > 0.");
  if (_Q <= 0.0)
    mooseError("InterstitialLoopGroupedNodalKernel requires Em_i_eV > 0.");
  if (_atomic_volume <= 0.0)
    mooseError("InterstitialLoopGroupedNodalKernel requires atomic_volume > 0.");
  if (_burgers_vector <= 0.0)
    mooseError(
        "InterstitialLoopGroupedNodalKernel requires lattice_parameter > 0 or burgers_vector > 0.");
  if (_lattice_parameter <= 0.0)
    mooseError("InterstitialLoopGroupedNodalKernel requires lattice_parameter > 0.");
  if (_formation_energy_eV <= 0.0)
    mooseError("InterstitialLoopGroupedNodalKernel requires formation_energy_eV > 0.");
  if (_energy_exponent <= 0.0)
    mooseError("InterstitialLoopGroupedNodalKernel requires loop_energy_exponent > 0.");
  if (_energy_coefficient_eV <= 0.0)
    mooseError("InterstitialLoopGroupedNodalKernel requires loop_energy_coefficient_eV > 0.");
  // The opposite-absorption parameters only need to be physically valid when the opposite
  // coupling is actually active: when opposite_coupled_var is not provided, this kernel reduces
  // exactly to its pre-existing self-reaction-only behavior (see validParams), and the unused
  // defaults (D0_v = 0, etc.) are never evaluated by any residual/Jacobian path -- all of the
  // Piece A/Piece B terms are gated on _opposite_coupling_active below.
  if (_opposite_coupling_active)
  {
    if (_opposite_D0 <= 0.0)
      mooseError("InterstitialLoopGroupedNodalKernel requires D0_v > 0 when opposite_coupled_var "
                 "is provided.");
    if (_opposite_Q <= 0.0)
      mooseError("InterstitialLoopGroupedNodalKernel requires Em_v_eV > 0 when "
                 "opposite_coupled_var is provided.");
    if (_opposite_group_nonnegative_tolerance_factor < 0.0)
      mooseError("InterstitialLoopGroupedNodalKernel requires "
                 "opposite_group_nonnegative_tolerance_factor >= 0.");
  }
}

template <bool is_ad>
Real
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::radius(unsigned int n) const
{
  return std::sqrt(static_cast<Real>(n) * _atomic_volume / (libMesh::pi * _burgers_vector));
}

template <bool is_ad>
Real
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::formationEnergy(unsigned int n) const
{
  return _energy_coefficient_eV * std::pow(static_cast<Real>(n), _energy_exponent) +
         _energy_constant_eV;
}

template <bool is_ad>
Real
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::loopBias(unsigned int n) const
{
  // Z_i(n) interpolates from the small-loop bias z_i to the large-loop bias z_loop_i as n grows,
  // following the dislocation-loop capture-efficiency model of Wiedersich/Jourdan-style AMSS
  // rate theory (DD_AMSS_Z_I, DD_AMSS_Z_LOOP_I, DD_AMSS_ALPHA_LOOP_I in the reference model).
  Real z = std::sqrt(_burgers_vector / (8.0 * libMesh::pi * _lattice_parameter)) * _z_loop_i - _z_i;
  z *= std::pow(1.0 / static_cast<Real>(n), 0.5 * _alpha_loop_i);
  return z + _z_i;
}

template <bool is_ad>
Real
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::monomerDiffusivity() const
{
  return _D0 * std::exp(-_Q / (loop_kB * _temperature));
}

template <bool is_ad>
Real
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::beta(unsigned int n) const
{
  if (n < _beta_cache.size())
    return _beta_cache[n];

  return 2.0 * libMesh::pi * radius(n) * loopBias(n) * monomerDiffusivity() / _atomic_volume;
}

template <bool is_ad>
Real
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::alpha(unsigned int n) const
{
  if (n < _alpha_cache.size())
    return _alpha_cache[n];

  if (n <= 1)
    mooseError("alpha(n) is only defined for loop sizes n >= 2.");

  const Real binding_energy_eV =
      _formation_energy_eV - (formationEnergy(n) - formationEnergy(n - 1));
  return beta(n - 1) * std::exp(-binding_energy_eV * loop_eV / (loop_kB * _temperature));
}

template <bool is_ad>
Real
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::oppositeMonomerDiffusivity() const
{
  return _opposite_D0 * std::exp(-_opposite_Q / (loop_kB * _temperature));
}

template <bool is_ad>
Real
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::biasVacancy(unsigned int n) const
{
  // Matches VacancyVoidGroupedNodalKernel::oppositeBias(), but evaluated with THIS kernel's own
  // burgers_vector/lattice_parameter (loopBiasVacancy() in the reference, same functional form as
  // loopBias() but with the vacancy-side bias parameters z_v/z_loop_v/alpha_loop_v).
  Real z = std::sqrt(_burgers_vector / (8.0 * libMesh::pi * _lattice_parameter)) * _z_loop_v - _z_v;
  z *= std::pow(1.0 / static_cast<Real>(n), 0.5 * _alpha_loop_v);
  return z + _z_v;
}

template <bool is_ad>
Real
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::betaOpposite(unsigned int n) const
{
  if (n < _beta_opposite_cache.size())
    return _beta_opposite_cache[n];

  return 2.0 * libMesh::pi * radius(n) * biasVacancy(n) * oppositeMonomerDiffusivity() /
         _atomic_volume;
}

template <bool is_ad>
Real
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::oppositeRadius(unsigned int n) const
{
  return std::cbrt(3.0 * static_cast<Real>(n) * _atomic_volume / (4.0 * libMesh::pi));
}

template <bool is_ad>
Real
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::oppositeBetaOpposite(unsigned int n) const
{
  if (n < _opposite_beta_opposite_cache.size())
    return _opposite_beta_opposite_cache[n];

  const Real r = oppositeRadius(n);
  return 4.0 * libMesh::pi * r * monomerDiffusivity() /
         ((1.0 + _lattice_parameter / r) * _atomic_volume);
}

template <bool is_ad>
void
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::ensureOppositeCoefficientCache() const
{
  if (_opposite_cache_size == _opposite_num_cluster_sizes)
    return;

  std::vector<Real> beta_opposite_values(_opposite_num_cluster_sizes + 1, 0.0);
  const Real diffusivity = monomerDiffusivity();
  for (unsigned int n = 1; n <= _opposite_num_cluster_sizes; ++n)
  {
    const Real r = oppositeRadius(n);
    beta_opposite_values[n] =
        4.0 * libMesh::pi * r * diffusivity / ((1.0 + _lattice_parameter / r) * _atomic_volume);
  }

  _opposite_beta_opposite_cache = std::move(beta_opposite_values);
  _opposite_cache_size = _opposite_num_cluster_sizes;
}

template <bool is_ad>
void
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::ensureCoefficientCache() const
{
  if (_cache_size == _num_cluster_sizes)
    return;

  // Compute into local vectors first: beta(n)/alpha(n)/betaOpposite(n) short-circuit to a cache
  // lookup once the caches have the right size, so filling them in place while they are already
  // sized would just read back the zeros they were assigned.
  std::vector<Real> beta_values(_num_cluster_sizes + 1, 0.0);
  std::vector<Real> alpha_values(_num_cluster_sizes + 1, 0.0);
  std::vector<Real> beta_opposite_values(_num_cluster_sizes + 1, 0.0);
  const Real diffusivity = monomerDiffusivity();
  const Real opposite_diffusivity = oppositeMonomerDiffusivity();

  for (unsigned int n = 1; n <= _num_cluster_sizes; ++n)
  {
    beta_values[n] = 2.0 * libMesh::pi * radius(n) * loopBias(n) * diffusivity / _atomic_volume;
    beta_opposite_values[n] =
        2.0 * libMesh::pi * radius(n) * biasVacancy(n) * opposite_diffusivity / _atomic_volume;
  }
  for (unsigned int n = 2; n <= _num_cluster_sizes; ++n)
  {
    const Real binding_energy_eV =
        _formation_energy_eV - (formationEnergy(n) - formationEnergy(n - 1));
    alpha_values[n] =
        beta_values[n - 1] * std::exp(-binding_energy_eV * loop_eV / (loop_kB * _temperature));
  }

  _beta_cache = std::move(beta_values);
  _alpha_cache = std::move(alpha_values);
  _beta_opposite_cache = std::move(beta_opposite_values);
  _cache_size = _num_cluster_sizes;
}

template <bool is_ad>
template <typename Scalar>
void
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::fillConcentrationCache(
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
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::fillFluxCache(
    const std::vector<Scalar> & concentration,
    const std::vector<Scalar> & opposite_concentration,
    Scalar opposite_monomer,
    std::vector<Scalar> & flux) const
{
  // combined_flux[n] = self_flux[n] - K[n+1], where K[n] = betaOpposite(n)*opposite_monomer*C_n is
  // the rate at which a vacancy monomer absorbed onto a loop of size n shrinks it to n-1
  // (oppositeShrinkFluxAt in the reference). See VacancyVoidGroupedNodalKernel::fillFluxCache for
  // the full derivation; this is the mirror image.
  flux.assign(_num_cluster_sizes + 1, Scalar(0.0));
  const Scalar c1 = concentration[1];
  for (unsigned int n = 1; n < _num_cluster_sizes; ++n)
  {
    const Scalar self_flux = beta(n) * c1 * concentration[n] - alpha(n + 1) * concentration[n + 1];
    const Scalar shrink_flux =
        n + 1 >= 2 ? betaOpposite(n + 1) * opposite_monomer * concentration[n + 1] : Scalar(0.0);
    flux[n] = self_flux - shrink_flux;
  }
}

template <bool is_ad>
template <typename Scalar>
void
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::computeGroupedResidualFromState(
    const std::vector<Scalar> & state,
    const std::vector<Scalar> & opposite_state,
    std::vector<Scalar> & residual_values) const
{
  std::vector<Scalar> concentration;
  fillConcentrationCache(state, concentration);

  std::vector<Scalar> opposite_concentration(_opposite_num_cluster_sizes + 1, Scalar(0.0));
  for (unsigned int n = 1; n <= _opposite_num_cluster_sizes; ++n)
    opposite_concentration[n] = _opposite_layout.concentration(
        opposite_state, n, _opposite_enforce_group_nonnegative,
        _opposite_group_nonnegative_tolerance_factor);
  const Scalar opposite_monomer = opposite_concentration[1];

  std::vector<Scalar> flux;
  fillFluxCache(concentration, opposite_concentration, opposite_monomer, flux);

  residual_values.assign(_layout.componentCount(), Scalar(0.0));

  // Monomer: unlike the explicit/grouped sizes, the monomer equation keeps the PLAIN self-flux
  // (not combined_flux) and gets one extra +K[2] term, plus the Piece B drain summing the entire
  // coupled void array. See VacancyVoidGroupedNodalKernel::computeGroupedResidualFromState for the
  // full derivation; this is the mirror image.
  Scalar self_monomer_flux = 2.0 * (beta(1) * concentration[1] * concentration[1] -
                                    alpha(2) * concentration[2]);
  for (unsigned int n = 2; n < _num_cluster_sizes; ++n)
    self_monomer_flux += beta(n) * concentration[1] * concentration[n] -
                         alpha(n + 1) * concentration[n + 1];
  const Scalar shrink_to_monomer =
      _num_cluster_sizes >= 2 ? betaOpposite(2) * opposite_monomer * concentration[2] : Scalar(0.0);

  Scalar opposite_monomer_loss = Scalar(0.0);
  for (unsigned int n = 2; n <= _opposite_num_cluster_sizes; ++n)
    opposite_monomer_loss +=
        oppositeBetaOpposite(n) * concentration[1] * opposite_concentration[n];

  residual_values[_layout.monomerComponent()] =
      -(_generation - _sink * concentration[1] - self_monomer_flux + shrink_to_monomer) +
      opposite_monomer_loss;

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
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::jacobianSetup()
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
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::computeQpResidual(
    GenericRealEigenVector<is_ad> & residual)
{
  const auto n_comp = _u[_qp].size();
  if (n_comp != _layout.componentCount())
    mooseError("InterstitialLoopGroupedNodalKernel expected ",
               _layout.componentCount(),
               " array components from the grouping layout, but the variable has ",
               n_comp,
               ".");

  ensureCoefficientCache();
  ensureOppositeCoefficientCache();
  std::vector<GenericReal<is_ad>> state(n_comp, GenericReal<is_ad>(0.0));
  for (unsigned int i = 0; i < n_comp; ++i)
    state[i] = _u[_qp](i);

  const auto opposite_n_comp = _opposite_array[_qp].size();
  if (opposite_n_comp != _opposite_layout.componentCount())
    mooseError("InterstitialLoopGroupedNodalKernel expected ",
               _opposite_layout.componentCount(),
               " array components from the opposite grouping layout, but opposite_coupled_var "
               "has ",
               opposite_n_comp,
               ".");
  std::vector<GenericReal<is_ad>> opposite_state(opposite_n_comp, GenericReal<is_ad>(0.0));
  for (unsigned int i = 0; i < opposite_n_comp; ++i)
    opposite_state[i] = _opposite_array[_qp](i);

  std::vector<GenericReal<is_ad>> residual_values;
  computeGroupedResidualFromState(state, opposite_state, residual_values);
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
GenericInterstitialLoopGroupedNodalKernelTempl<false>::computeQpJacobian()
{
  const auto n_comp = static_cast<unsigned int>(_u[_qp].size());
  if (n_comp != _layout.componentCount())
    mooseError("InterstitialLoopGroupedNodalKernel expected ",
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

  ensureOppositeCoefficientCache();
  const Real opposite_monomer = _opposite_array[_qp](0);
  Real opposite_coefficient_sum = 0.0;
  for (unsigned int n = 2; n <= _opposite_num_cluster_sizes; ++n)
  {
    const Real opposite_c_n = _opposite_layout.concentration(
        _opposite_array[_qp], n, _opposite_enforce_group_nonnegative,
        _opposite_group_nonnegative_tolerance_factor);
    opposite_coefficient_sum += oppositeBetaOpposite(n) * opposite_c_n;
  }

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
  // Plain self-flux derivative (alpha only): used ONLY for the monomer row's self_monomer_flux
  // term, which by construction (see computeGroupedResidualFromState) does not include the
  // opposite-absorption K terms -- those are added to the monomer row as separate explicit terms
  // below, not folded into this helper.
  auto addFluxDerivative = [&](unsigned int n, Real weight)
  {
    // J_0 and J_N vanish at the ends of the truncated size space
    if (n == 0 || n >= _num_cluster_sizes)
      return;
    addConcentrationDerivative(1, weight * beta(n) * concentration[n]);
    addConcentrationDerivative(n, weight * beta(n) * c1);
    addConcentrationDerivative(n + 1, -weight * alpha(n + 1));
  };
  // Combined-flux derivative (alpha + opposite-shrink): used for explicit-size and grouped-bin
  // rows, where the residual's flux[n] means combined_flux[n] = self_flux[n] -
  // betaOpposite(n+1)*opposite_monomer*C_{n+1}. See VacancyVoidGroupedNodalKernel for the full
  // reasoning; this is the mirror image.
  auto addCombinedFluxDerivative = [&](unsigned int n, Real weight)
  {
    if (n == 0 || n >= _num_cluster_sizes)
      return;
    addConcentrationDerivative(1, weight * beta(n) * concentration[n]);
    addConcentrationDerivative(n, weight * beta(n) * c1);
    const Real effective_alpha = alpha(n + 1) + betaOpposite(n + 1) * opposite_monomer;
    addConcentrationDerivative(n + 1, -weight * effective_alpha);
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

  // Monomer: F = -(G_1 - k_s*C_1 - self_monomer_flux + shrink_to_monomer) + opposite_monomer_loss
  addEntry(_layout.monomerComponent(), _sink);
  for (unsigned int n = 1; n < _num_cluster_sizes; ++n)
    addFluxDerivative(n, n == 1 ? 2.0 : 1.0);
  // -shrink_to_monomer = -betaOpposite(2)*opposite_monomer*C_2: intra-variable derivative wrt C_2
  if (_num_cluster_sizes >= 2)
    addConcentrationDerivative(2, -betaOpposite(2) * opposite_monomer);
  // +opposite_monomer_loss = C_1 * sum_n oppositeBetaOpposite(n)*C_opposite,n: intra-variable
  // derivative wrt C_1 is the (fixed, within this call) coefficient sum itself.
  addConcentrationDerivative(1, opposite_coefficient_sum);
  setRow(_layout.monomerComponent());

  // Explicit sizes: F = -(J_{n-1} - J_n), J now meaning combined_flux
  for (unsigned int n = 2; n <= _layout.explicitMax(); ++n)
  {
    addCombinedFluxDerivative(n - 1, -1.0);
    addCombinedFluxDerivative(n, 1.0);
    setRow(_layout.explicitComponent(n));
  }

  for (const auto & bin : _layout.groups())
  {
    // L0: F = -(J_left - J_right) / width
    addCombinedFluxDerivative(bin.start - 1, -1.0 / bin.width);
    addCombinedFluxDerivative(bin.end, 1.0 / bin.width);
    setRow(bin.l0_component);

    // L1: F = (J_left - 2*J_avg + J_right) * (width - 1) / (2*sigma2*width)
    if (bin.sigma2 > 0.0 && bin.width > 1.0)
    {
      const Real factor = (bin.width - 1.0) / (2.0 * bin.sigma2 * bin.width);
      addCombinedFluxDerivative(bin.start - 1, factor);
      for (unsigned int n = bin.start; n < bin.end; ++n)
        addCombinedFluxDerivative(n, -2.0 * factor / (bin.width - 1.0));
      addCombinedFluxDerivative(bin.end, factor);
      setRow(bin.l1_component);
    }
  }
}

template <>
void
GenericInterstitialLoopGroupedNodalKernelTempl<true>::computeQpJacobian()
{
  mooseError("Internal error: computeQpJacobian should never be called for the AD version");
}

// Cross-variable (opposite-array) Jacobian for the opposite-absorption terms. Non-AD only, same
// reasoning as PointDefectRecombinationNodalKernel/VacancyVoidGroupedNodalKernel: the AD variant's
// computeOffDiagJacobian only re-runs the AD residual graph, which never touched the opposite
// array's own AD derivatives (coupledArrayValue strips them), so it would silently assemble a zero
// block here.
template <>
void
GenericInterstitialLoopGroupedNodalKernelTempl<false>::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar != _opposite_var_num)
    return;

  ensureCoefficientCache();
  ensureOppositeCoefficientCache();

  const RealEigenVector & state = _u[_qp];
  std::vector<Real> concentration(_num_cluster_sizes + 1, 0.0);
  for (unsigned int n = 1; n <= _num_cluster_sizes; ++n)
    concentration[n] = _layout.concentration(
        state, n, _enforce_group_nonnegative, _group_nonnegative_tolerance_factor);

  const auto opposite_n_comp = static_cast<unsigned int>(_opposite_array[_qp].size());
  std::vector<Real> opposite_a0(_opposite_num_cluster_sizes + 1, 1.0);
  std::vector<Real> opposite_a1(_opposite_num_cluster_sizes + 1, 0.0);
  std::vector<Real> opposite_concentration(_opposite_num_cluster_sizes + 1, 0.0);
  for (unsigned int n = 1; n <= _opposite_num_cluster_sizes; ++n)
  {
    if (_opposite_layout.isExplicitSize(n))
    {
      opposite_concentration[n] = _opposite_array[_qp](_opposite_layout.explicitComponent(n));
      continue;
    }
    const auto & bin = _opposite_layout.groupForSize(n);
    const Real l0 = _opposite_array[_qp](bin.l0_component);
    const Real l1 = _opposite_array[_qp](bin.l1_component);
    ClusterGroupingLayout::reconstructionCoefficients(bin,
                                                      n,
                                                      l0,
                                                      l1,
                                                      _opposite_enforce_group_nonnegative,
                                                      _opposite_group_nonnegative_tolerance_factor,
                                                      opposite_a0[n],
                                                      opposite_a1[n]);
    opposite_concentration[n] = opposite_a0[n] * l0 + opposite_a1[n] * l1;
  }

  std::vector<Real> row_values(opposite_n_comp, 0.0);
  std::vector<bool> row_touched(opposite_n_comp, false);
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
  // d(combined_flux[n])/d(opposite monomer) = -betaOpposite(n+1)*C_{n+1}(own): a single entry at
  // the opposite array's monomer component (0).
  auto addFluxOffDiagDerivative = [&](unsigned int n, Real weight)
  {
    if (n == 0 || n >= _num_cluster_sizes)
      return;
    addEntry(0, -weight * betaOpposite(n + 1) * concentration[n + 1]);
  };
  auto addOppositeConcentrationDerivative = [&](unsigned int n, Real value)
  {
    if (_opposite_layout.isExplicitSize(n))
      addEntry(_opposite_layout.explicitComponent(n), value);
    else
    {
      const auto & bin = _opposite_layout.groupForSize(n);
      addEntry(bin.l0_component, value * opposite_a0[n]);
      addEntry(bin.l1_component, value * opposite_a1[n]);
    }
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

  // Monomer row: d(-shrink_to_monomer)/d(opposite monomer) = -betaOpposite(2)*C_2(own); plus
  // d(opposite_monomer_loss)/d(opposite C_n) = C_1(own) * oppositeBetaOpposite(n) for every
  // opposite size n >= 2 (dense row against the opposite array, mirroring the existing
  // intra-variable Jacobian's own dense monomer row against its own array).
  if (_num_cluster_sizes >= 2)
    addEntry(0, -betaOpposite(2) * concentration[2]);
  for (unsigned int n = 2; n <= _opposite_num_cluster_sizes; ++n)
    addOppositeConcentrationDerivative(n, concentration[1] * oppositeBetaOpposite(n));
  setRow(_layout.monomerComponent());

  // Explicit sizes and grouped bins: only the single opposite-monomer entry from
  // addFluxOffDiagDerivative; no other dependence on the opposite array.
  for (unsigned int n = 2; n <= _layout.explicitMax(); ++n)
  {
    addFluxOffDiagDerivative(n - 1, -1.0);
    addFluxOffDiagDerivative(n, 1.0);
    setRow(_layout.explicitComponent(n));
  }
  for (const auto & bin : _layout.groups())
  {
    addFluxOffDiagDerivative(bin.start - 1, -1.0 / bin.width);
    addFluxOffDiagDerivative(bin.end, 1.0 / bin.width);
    setRow(bin.l0_component);

    if (bin.sigma2 > 0.0 && bin.width > 1.0)
    {
      const Real factor = (bin.width - 1.0) / (2.0 * bin.sigma2 * bin.width);
      addFluxOffDiagDerivative(bin.start - 1, factor);
      for (unsigned int n = bin.start; n < bin.end; ++n)
        addFluxOffDiagDerivative(n, -2.0 * factor / (bin.width - 1.0));
      addFluxOffDiagDerivative(bin.end, factor);
      setRow(bin.l1_component);
    }
  }
}

template <>
void
GenericInterstitialLoopGroupedNodalKernelTempl<true>::computeQpOffDiagJacobian(unsigned int)
{
  // Intentionally a no-op: see the class comment and PointDefectRecombinationNodalKernel for why
  // the AD variant cannot assemble this cross-variable block.
}

template class GenericInterstitialLoopGroupedNodalKernelTempl<false>;
template class GenericInterstitialLoopGroupedNodalKernelTempl<true>;
