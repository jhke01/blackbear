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

#include "GenericClusterDynamicsGroupedNodalKernel.h"
#include "ClusterDynamicsSolverCheck.h"
#include "FEProblemBase.h"
#include "NonlinearSystemBase.h"

#include "libmesh/libmesh_common.h"

#include <cstddef>

namespace
{
constexpr Real grouped_kB = 1.380649e-23;
constexpr Real grouped_eV = 1.602176634e-19;

std::vector<Real>
convertGroupedEvVectorToJ(const std::vector<Real> & values)
{
  std::vector<Real> converted(values.size());
  for (std::size_t i = 0; i < values.size(); ++i)
    converted[i] = values[i] * grouped_eV;
  return converted;
}

// 4*exp(-E_f/kT), matching DefectDynamicsAction::monomerConcentration()'s identical formula and
// physical constants exactly. Guarded against temperature == 0 (rate_model = simple, where this
// value is never used) so the member initializer never divides by zero.
Real
vacancyThermalConcentrationGrouped(Real formation_energy_eV, Real temperature)
{
  if (temperature <= 0.0)
    return 0.0;
  return 4.0 * std::exp(-formation_energy_eV * grouped_eV / (grouped_kB * temperature));
}
} // namespace

registerMooseObject("BlackBearApp", ClusterDynamicsGroupedNodalKernel);
registerMooseObject("BlackBearApp", ADClusterDynamicsGroupedNodalKernel);

template <bool is_ad>
InputParameters
GenericClusterDynamicsGroupedNodalKernelTempl<is_ad>::validParams()
{
  InputParameters params = GenericArrayNodalKernel<is_ad>::validParams();
  MooseEnum rate_model("simple interfacial_energy", "simple");
  MooseEnum diffusivity_model("constant arrhenius", "constant");
  MooseEnum binding_energy_model("interfacial_energy capillary binding_energy_table "
                                 "formation_energy_table",
                                 "interfacial_energy");
  MooseEnum group_spacing("linear geometric", "linear");
  params.addClassDescription(
      "Grouped cluster dynamics nodal kernel using explicit small sizes and an "
      "L0/L1 grouped tail for large clusters.");
  params.addParam<Real>("generation",
                        0.0,
                        "Monomer generation rate G_1. Applies to irradiation-produced point "
                        "defects; leave at zero for solute clustering such as Cu precipitation.");
  params.addParam<Real>("sink",
                        0.0,
                        "Linear sink coefficient k_s for monomer loss. Applies to point defects; "
                        "leave at zero for solute clustering such as Cu precipitation.");
  params.addRequiredParam<unsigned int>(
      "num_cluster_sizes",
      "Maximum physical cluster size represented by the grouped distribution.");
  params.addParam<bool>(
      "use_grouping", true, "Whether to activate grouped-tail storage above explicit_max.");
  params.addParam<unsigned int>(
      "explicit_max", 100, "Maximum cluster size retained explicitly before grouped bins.");
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
  params.addParam<MooseEnum>("rate_model",
                             rate_model,
                             "Rate-coefficient model. 'simple' uses beta_n = beta0*n^(1/3) and "
                             "alpha_n = alpha0*n^(1/3); 'interfacial_energy' derives coefficients "
                             "from cluster geometry, interfacial energy, and detailed balance.");
  params.addParam<MooseEnum>(
      "diffusivity_model",
      diffusivity_model,
      "Diffusivity model used only when rate_model = interfacial_energy. "
      "'constant' uses monomer_diffusivity directly, while 'arrhenius' computes "
      "the diffusivity from D0, Q_eV, and temperature.");
  params.addParam<MooseEnum>(
      "binding_energy_model",
      binding_energy_model,
      "Binding-energy model used by detailed-balance emission when rate_model = "
      "interfacial_energy. 'interfacial_energy' uses Omega, DeltaS, sigma, and atomic "
      "volume; 'capillary' uses A - B*(n^p - (n-1)^p), where A may be fixed or "
      "temperature dependent and B may be supplied directly or computed from sigma and "
      "atomic_volume; 'binding_energy_table' reads direct binding energies for n >= 2; "
      "and 'formation_energy_table' computes E_b(n) = E_f(n-1) + E_f(1) - E_f(n).");
  params.addParam<Real>("beta0", 0.0, "Base absorption coefficient for rate_model = simple.");
  params.addParam<Real>("alpha0", 0.0, "Base emission coefficient for rate_model = simple.");
  params.addParam<Real>(
      "temperature", 0.0, "Temperature T [K] for the interfacial-energy rate model.");
  params.addParam<Real>("monomer_diffusivity",
                        0.0,
                        "Monomer diffusion coefficient [m^2/s] for diffusivity_model = constant.");
  params.addParam<Real>(
      "D0", 0.0, "Diffusion prefactor D0 [m^2/s] for diffusivity_model = arrhenius.");
  params.addParam<Real>("Q_eV", 0.0, "Activation energy Q [eV] for diffusivity_model = arrhenius.");
  params.addParam<Real>("radiation_enhanced_factor",
                        1.0,
                        "Multiplicative radiation-enhanced factor applied to the monomer "
                        "diffusivity for the interfacial-energy model.");
  params.addCoupledVar(
      "vacancy_coupled_var",
      0.0,
      "Optional DefectDynamics vacancy grouped array variable. When supplied, the monomer "
      "diffusivity is additionally multiplied by the dynamic vacancy supersaturation ratio "
      "C_v(1,t)/C_v_thermal, in place of (or on top of) radiation_enhanced_factor. "
      "C_v_thermal is computed from vacancy_formation_energy_eV at this kernel's own "
      "temperature.");
  params.addParam<Real>(
      "vacancy_formation_energy_eV",
      0.0,
      "Vacancy monomer formation energy E_f [eV], required when vacancy_coupled_var is given -- "
      "used to compute C_v_thermal = 4*exp(-E_f/kT) at this kernel's own temperature, matching "
      "DefectDynamicsAction::monomerConcentration()'s identical formula and physical constants.");
  params.addParam<Real>(
      "sigma",
      0.0,
      "Interfacial energy sigma [J/m^2] for binding_energy_model = interfacial_energy "
      "or capillary when the capillary surface coefficient is computed from sigma.");
  params.addParam<Real>(
      "atomic_volume", 0.0, "Atomic volume V_at [m^3] for the interfacial-energy model.");
  params.addParam<Real>(
      "Omega_kB_K", 0.0, "Enthalpy term Omega [k_B*K] for the interfacial-energy model.");
  params.addParam<Real>(
      "DeltaS_kB",
      0.0,
      "Non-configurational entropy term DeltaS [k_B] for the interfacial-energy model.");
  params.addParam<Real>("binding_energy_constant_eV",
                        0.0,
                        "Constant A [eV] in binding_energy_model = capillary: "
                        "E_b(n) = A - B*(n^p - (n-1)^p).");
  params.addParam<Real>(
      "binding_energy_enthalpy_kB_K",
      0.0,
      "Enthalpy-like term H [k_B*K] for the temperature-dependent constant in "
      "binding_energy_model = capillary: E_b(n) = (H - T*S)*k_B - B*(n^p - (n-1)^p).");
  params.addParam<Real>("binding_energy_entropy_kB",
                        0.0,
                        "Entropy-like term S [k_B] for the temperature-dependent constant in "
                        "binding_energy_model = capillary.");
  params.addParam<Real>(
      "binding_energy_coefficient_eV",
      0.0,
      "Direct coefficient B [eV] in binding_energy_model = capillary. Leave this at zero "
      "and provide sigma to compute B from interfacial energy instead.");
  params.addParam<Real>(
      "binding_energy_exponent", 2.0 / 3.0, "Exponent p in binding_energy_model = capillary.");
  params.addParam<std::vector<Real>>(
      "binding_energy_table_eV",
      std::vector<Real>(),
      "Direct binding energies [eV] for binding_energy_model = binding_energy_table. "
      "Entry 0 corresponds to cluster size n = 2.");
  params.addParam<std::vector<Real>>(
      "formation_energy_table_eV",
      std::vector<Real>(),
      "Formation energies [eV] for binding_energy_model = formation_energy_table. "
      "Entry 0 corresponds to cluster size n = 1.");
  return params;
}

template <bool is_ad>
GenericClusterDynamicsGroupedNodalKernelTempl<is_ad>::GenericClusterDynamicsGroupedNodalKernelTempl(
    const InputParameters & parameters)
  : GenericArrayNodalKernel<is_ad>(parameters),
    _generation(this->template getParam<Real>("generation")),
    _sink(this->template getParam<Real>("sink")),
    _rate_model(this->template getParam<MooseEnum>("rate_model").template getEnum<RateModel>()),
    _beta0(this->template getParam<Real>("beta0")),
    _alpha0(this->template getParam<Real>("alpha0")),
    _temperature(this->template getParam<Real>("temperature")),
    _monomer_diffusivity(this->template getParam<Real>("monomer_diffusivity")),
    _diffusivity_model(this->template getParam<MooseEnum>("diffusivity_model")
                           .template getEnum<DiffusivityModel>()),
    _binding_energy_model(this->template getParam<MooseEnum>("binding_energy_model")
                              .template getEnum<BindingEnergyModel>()),
    _D0(this->template getParam<Real>("D0")),
    _Q(this->template getParam<Real>("Q_eV") * grouped_eV),
    _radiation_enhanced_factor(this->template getParam<Real>("radiation_enhanced_factor")),
    _vacancy_array(this->coupledArrayValue("vacancy_coupled_var")),
    _vacancy_coupling_active(this->isCoupled("vacancy_coupled_var")),
    _vacancy_var_num(this->coupled("vacancy_coupled_var")),
    _vacancy_c_thermal(vacancyThermalConcentrationGrouped(
        this->template getParam<Real>("vacancy_formation_energy_eV"),
        this->template getParam<Real>("temperature"))),
    _sigma(this->template getParam<Real>("sigma")),
    _atomic_volume(this->template getParam<Real>("atomic_volume")),
    _Omega(this->template getParam<Real>("Omega_kB_K") * grouped_kB),
    _DeltaS(this->template getParam<Real>("DeltaS_kB") * grouped_kB),
    _binding_energy_constant(this->template getParam<Real>("binding_energy_constant_eV") *
                             grouped_eV),
    _binding_energy_enthalpy(this->template getParam<Real>("binding_energy_enthalpy_kB_K") *
                             grouped_kB),
    _binding_energy_entropy(this->template getParam<Real>("binding_energy_entropy_kB") *
                            grouped_kB),
    _binding_energy_coefficient(this->template getParam<Real>("binding_energy_coefficient_eV") *
                                grouped_eV),
    _binding_energy_exponent(this->template getParam<Real>("binding_energy_exponent")),
    _binding_energy_table(convertGroupedEvVectorToJ(
        this->template getParam<std::vector<Real>>("binding_energy_table_eV"))),
    _formation_energy_table(convertGroupedEvVectorToJ(
        this->template getParam<std::vector<Real>>("formation_energy_table_eV"))),
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
  if (_vacancy_coupling_active &&
      this->template getParam<Real>("vacancy_formation_energy_eV") <= 0.0)
    mooseError("ClusterDynamicsGroupedNodalKernel requires vacancy_formation_energy_eV > 0 when "
              "vacancy_coupled_var is given.");

  // radiation_enhanced_factor and vacancy_coupled_var both multiply the SAME monomer
  // diffusivity: a nonzero vacancy supersaturation ratio (S_v = C_v(1,t)/C_v_thermal) is itself
  // often many orders of magnitude under irradiation, so leaving radiation_enhanced_factor at a
  // value calibrated to approximate that same enhancement statically (its usual role before this
  // coupling existed) double-counts the effect rather than combining two independent physical
  // mechanisms.
  if (_vacancy_coupling_active && this->parameters().isParamSetByUser("radiation_enhanced_factor") &&
      _radiation_enhanced_factor != 1.0)
    this->paramWarning(
        "radiation_enhanced_factor",
        "radiation_enhanced_factor and vacancy_coupled_var both scale the same monomer "
        "diffusivity; leaving radiation_enhanced_factor at a value other than 1.0 while "
        "vacancy_coupled_var is active double-counts the enhancement instead of replacing a "
        "static placeholder with the dynamic vacancy supersaturation ratio. Set "
        "radiation_enhanced_factor = 1.0 (its default) unless the two are deliberately meant to "
        "compound.");

  // Grouped and explicit components are coupled through the monomer and the inter-bin fluxes,
  // so the Jacobian contains entries outside the default diagonal-block sparsity pattern.
  if (!this->_fe_problem.useHashTableMatrixAssembly())
    mooseError(
        "ClusterDynamicsGroupedNodalKernel requires Problem/use_hash_table_matrix_assembly = true");
  if (this->_tid == 0)
    ClusterDynamics::checkSolverSetup(*this,
                                      this->_fe_problem,
                                      this->_var.sys().number(),
                                      this->_var.number(),
                                      this->_var.count());

  if (_group_nonnegative_tolerance_factor < 0.0)
    mooseError(
        "ClusterDynamicsGroupedNodalKernel requires group_nonnegative_tolerance_factor >= 0.");

  if (_rate_model == RateModel::SIMPLE)
  {
    if (_beta0 <= 0.0)
      mooseError("ClusterDynamicsGroupedNodalKernel with rate_model = simple requires beta0 > 0.");
    if (_alpha0 <= 0.0)
      mooseError("ClusterDynamicsGroupedNodalKernel with rate_model = simple requires alpha0 > 0.");
  }
  else
  {
    if (_temperature <= 0.0)
      mooseError("ClusterDynamicsGroupedNodalKernel with rate_model = interfacial_energy requires "
                 "temperature > 0.");
    if (_diffusivity_model == DiffusivityModel::CONSTANT)
    {
      if (_monomer_diffusivity <= 0.0)
        mooseError("ClusterDynamicsGroupedNodalKernel with rate_model = interfacial_energy and "
                   "diffusivity_model = constant requires monomer_diffusivity > 0.");
    }
    else
    {
      if (_D0 <= 0.0)
        mooseError("ClusterDynamicsGroupedNodalKernel with rate_model = interfacial_energy and "
                   "diffusivity_model = arrhenius requires D0 > 0.");
      if (_Q <= 0.0)
        mooseError("ClusterDynamicsGroupedNodalKernel with rate_model = interfacial_energy and "
                   "diffusivity_model = arrhenius requires Q_eV > 0.");
    }
    if (_atomic_volume <= 0.0)
      mooseError("ClusterDynamicsGroupedNodalKernel with rate_model = interfacial_energy requires "
                 "atomic_volume > 0.");
    if (_radiation_enhanced_factor <= 0.0)
      mooseError("ClusterDynamicsGroupedNodalKernel with rate_model = interfacial_energy requires "
                 "radiation_enhanced_factor > 0.");
    if (_binding_energy_model == BindingEnergyModel::INTERFACIAL_ENERGY)
    {
      if (_sigma <= 0.0)
        mooseError("ClusterDynamicsGroupedNodalKernel with binding_energy_model = "
                   "interfacial_energy requires sigma > 0.");
      if (_Omega <= 0.0)
        mooseError("ClusterDynamicsGroupedNodalKernel with binding_energy_model = "
                   "interfacial_energy requires Omega_kB_K > 0.");
      if (_DeltaS <= 0.0)
        mooseError("ClusterDynamicsGroupedNodalKernel with binding_energy_model = "
                   "interfacial_energy requires DeltaS_kB > 0.");
    }
    else if (_binding_energy_model == BindingEnergyModel::CAPILLARY)
    {
      const bool has_fixed_constant = _binding_energy_constant != 0.0;
      const bool has_thermal_constant =
          _binding_energy_enthalpy != 0.0 || _binding_energy_entropy != 0.0;
      const bool has_direct_coefficient = _binding_energy_coefficient != 0.0;
      const bool has_surface_coefficient = _sigma != 0.0;

      if (has_fixed_constant && has_thermal_constant)
        mooseError("ClusterDynamicsGroupedNodalKernel with binding_energy_model = capillary "
                   "cannot use both binding_energy_constant_eV and binding_energy_enthalpy_kB_K/"
                   "binding_energy_entropy_kB.");
      if (has_direct_coefficient && has_surface_coefficient)
        mooseError("ClusterDynamicsGroupedNodalKernel with binding_energy_model = capillary "
                   "cannot use both binding_energy_coefficient_eV and sigma.");
      if (!has_fixed_constant && !has_thermal_constant && !has_direct_coefficient &&
          !has_surface_coefficient)
        mooseError("ClusterDynamicsGroupedNodalKernel with binding_energy_model = capillary "
                   "requires a constant term, a thermal H/S term, a direct coefficient, or "
                   "sigma.");
      if (has_surface_coefficient && _sigma <= 0.0)
        mooseError("ClusterDynamicsGroupedNodalKernel with binding_energy_model = capillary "
                   "requires sigma > 0 when sigma is used for the size coefficient.");
      if (_binding_energy_exponent <= 0.0)
        mooseError("ClusterDynamicsGroupedNodalKernel requires binding_energy_exponent > 0.");
    }
    else if (_binding_energy_model == BindingEnergyModel::BINDING_ENERGY_TABLE)
    {
      if (_binding_energy_table.empty())
        mooseError("ClusterDynamicsGroupedNodalKernel with binding_energy_model = "
                   "binding_energy_table requires binding_energy_table_eV.");
    }
    else if (_binding_energy_model == BindingEnergyModel::FORMATION_ENERGY_TABLE)
    {
      if (_formation_energy_table.size() < 2)
        mooseError("ClusterDynamicsGroupedNodalKernel with binding_energy_model = "
                   "formation_energy_table requires at least E_f(1) and E_f(2).");
    }
  }
}

template <bool is_ad>
Real
GenericClusterDynamicsGroupedNodalKernelTempl<is_ad>::beta(unsigned int n) const
{
  if (n < _beta_cache.size())
    // In interfacial_energy mode, _beta_cache holds beta(n)/monomerDiffusivity() (see
    // ensureCoefficientCache()): a pure function of n that never changes, so the live, possibly
    // solution-dependent diffusivity is applied here as a cheap O(1) multiply rather than baked
    // into an O(_num_cluster_sizes) cache rebuild on every call.
    return _rate_model == RateModel::SIMPLE ? _beta_cache[n] : _beta_cache[n] * monomerDiffusivity();

  if (_rate_model == RateModel::SIMPLE)
    return _beta0 * std::cbrt(static_cast<Real>(n));

  return 4.0 * libMesh::pi * (radius(1) + radius(n)) * monomerDiffusivity() / atomicVolume();
}

template <bool is_ad>
Real
GenericClusterDynamicsGroupedNodalKernelTempl<is_ad>::alpha(unsigned int n) const
{
  if (n < _alpha_cache.size())
    // Same shape-cache convention as beta() above: _alpha_cache[n] = alpha(n)/monomerDiffusivity()
    // in interfacial_energy mode.
    return _rate_model == RateModel::SIMPLE ? _alpha_cache[n] : _alpha_cache[n] * monomerDiffusivity();

  if (_rate_model == RateModel::SIMPLE)
    return _alpha0 * std::cbrt(static_cast<Real>(n));

  if (n <= 1)
    mooseError("alpha(n) is only defined for cluster sizes n >= 2 in interfacial_energy mode.");

  return beta(n - 1) * std::exp(-bindingEnergy(n) / (grouped_kB * _temperature));
}

template <bool is_ad>
Real
GenericClusterDynamicsGroupedNodalKernelTempl<is_ad>::atomicVolume() const
{
  return _atomic_volume;
}

template <bool is_ad>
Real
GenericClusterDynamicsGroupedNodalKernelTempl<is_ad>::baseDiffusivity() const
{
  if (_diffusivity_model == DiffusivityModel::CONSTANT)
    return _monomer_diffusivity;
  return _D0 * std::exp(-_Q / (grouped_kB * _temperature));
}

template <bool is_ad>
Real
GenericClusterDynamicsGroupedNodalKernelTempl<is_ad>::monomerDiffusivity() const
{
  Real diffusivity = _radiation_enhanced_factor * baseDiffusivity();
  if (_vacancy_coupling_active)
    diffusivity *= _vacancy_array[_qp](0) / _vacancy_c_thermal;
  return diffusivity;
}

template <bool is_ad>
Real
GenericClusterDynamicsGroupedNodalKernelTempl<is_ad>::monomerDiffusivityCouplingDerivative() const
{
  if (!_vacancy_coupling_active)
    return 0.0;
  return _radiation_enhanced_factor * baseDiffusivity() / _vacancy_c_thermal;
}

template <bool is_ad>
Real
GenericClusterDynamicsGroupedNodalKernelTempl<is_ad>::betaPrime(const unsigned int n) const
{
  if (!_vacancy_coupling_active || _rate_model == RateModel::SIMPLE)
    return 0.0;
  return 4.0 * libMesh::pi * (radius(1) + radius(n)) * monomerDiffusivityCouplingDerivative() /
         atomicVolume();
}

template <bool is_ad>
Real
GenericClusterDynamicsGroupedNodalKernelTempl<is_ad>::alphaPrime(const unsigned int n) const
{
  if (!_vacancy_coupling_active || _rate_model == RateModel::SIMPLE)
    return 0.0;
  if (n <= 1)
    mooseError(
        "alphaPrime(n) is only defined for cluster sizes n >= 2 in interfacial_energy mode.");
  // alpha(n) = beta(n-1)*exp(-E_b(n)/kT); E_b(n) does not depend on the coupled vacancy
  // concentration, so this is the same chain rule with beta replaced by betaPrime.
  return betaPrime(n - 1) * std::exp(-bindingEnergy(n) / (grouped_kB * _temperature));
}

template <bool is_ad>
Real
GenericClusterDynamicsGroupedNodalKernelTempl<is_ad>::radius(unsigned int n) const
{
  return std::cbrt(3.0 * static_cast<Real>(n) * atomicVolume() / (4.0 * libMesh::pi));
}

template <bool is_ad>
Real
GenericClusterDynamicsGroupedNodalKernelTempl<is_ad>::bindingEnergy(unsigned int n) const
{
  if (n <= 1)
    mooseError("bindingEnergy(n) is only defined for cluster sizes n >= 2.");

  if (_binding_energy_model == BindingEnergyModel::CAPILLARY)
  {
    const bool has_fixed_constant = _binding_energy_constant != 0.0;
    const bool has_direct_coefficient = _binding_energy_coefficient != 0.0;
    const Real constant_term =
        has_fixed_constant ? _binding_energy_constant
                           : _binding_energy_enthalpy - _temperature * _binding_energy_entropy;
    const Real size_coefficient =
        has_direct_coefficient
            ? _binding_energy_coefficient
            : std::cbrt(36.0 * libMesh::pi) * std::pow(atomicVolume(), 2.0 / 3.0) * _sigma;
    return constant_term -
           size_coefficient * (std::pow(static_cast<Real>(n), _binding_energy_exponent) -
                               std::pow(static_cast<Real>(n - 1), _binding_energy_exponent));
  }

  if (_binding_energy_model == BindingEnergyModel::BINDING_ENERGY_TABLE)
  {
    const auto index = n - 2;
    if (index >= _binding_energy_table.size())
      mooseError("binding_energy_table_eV does not contain a value for cluster size n = ", n, ".");
    return _binding_energy_table[index];
  }

  if (_binding_energy_model == BindingEnergyModel::FORMATION_ENERGY_TABLE)
  {
    if (n > _formation_energy_table.size())
      mooseError(
          "formation_energy_table_eV does not contain a value for cluster size n = ", n, ".");
    return _formation_energy_table[n - 2] + _formation_energy_table[0] -
           _formation_energy_table[n - 1];
  }

  const Real surface_prefactor =
      std::cbrt(36.0 * libMesh::pi) * std::pow(atomicVolume(), 2.0 / 3.0) * _sigma;
  return _Omega - _temperature * _DeltaS -
         surface_prefactor * (std::pow(static_cast<Real>(n), 2.0 / 3.0) -
                              std::pow(static_cast<Real>(n - 1), 2.0 / 3.0));
}

template <bool is_ad>
void
GenericClusterDynamicsGroupedNodalKernelTempl<is_ad>::ensureCoefficientCache() const
{
  // In interfacial_energy mode, beta(n) = shape(n) * monomerDiffusivity() and
  // alpha(n) = shape(n-1) * monomerDiffusivity() * exp(-Eb(n)/kT): the diffusivity factors out of
  // every cached value as one common multiplier. _beta_cache/_alpha_cache below therefore store
  // beta(n)/monomerDiffusivity() and alpha(n)/monomerDiffusivity() (that pure, n-only shape, never
  // the live diffusivity itself), so this cache depends only on _num_cluster_sizes and is rebuilt
  // exactly as rarely as before vacancy_coupled_var existed -- beta()/alpha() apply the current,
  // possibly solution-dependent monomerDiffusivity() as a cheap O(1) multiply at query time
  // instead.
  if (_cache_size == _num_cluster_sizes)
    return;

  _beta_cache.assign(_num_cluster_sizes + 1, 0.0);
  _alpha_cache.assign(_num_cluster_sizes + 1, 0.0);

  if (_rate_model == RateModel::SIMPLE)
  {
    for (unsigned int n = 1; n <= _num_cluster_sizes; ++n)
    {
      const Real n_to_one_third = std::cbrt(static_cast<Real>(n));
      _beta_cache[n] = _beta0 * n_to_one_third;
      if (n >= 2)
        _alpha_cache[n] = _alpha0 * n_to_one_third;
    }
  }
  else
  {
    const Real atomic_volume = atomicVolume();
    const Real inv_atomic_volume = 1.0 / atomic_volume;
    const Real radius_prefactor = std::cbrt(3.0 * atomic_volume / (4.0 * libMesh::pi));
    const Real r1 = radius_prefactor;
    const Real thermal_factor = grouped_kB * _temperature;

    std::vector<Real> radius_cache(_num_cluster_sizes + 1, 0.0);
    for (unsigned int n = 1; n <= _num_cluster_sizes; ++n)
    {
      radius_cache[n] = radius_prefactor * std::cbrt(static_cast<Real>(n));
      _beta_cache[n] = 4.0 * libMesh::pi * (r1 + radius_cache[n]) * inv_atomic_volume;
    }

    for (unsigned int n = 2; n <= _num_cluster_sizes; ++n)
      _alpha_cache[n] = _beta_cache[n - 1] * std::exp(-bindingEnergy(n) / thermal_factor);
  }

  _cache_size = _num_cluster_sizes;
}

template <bool is_ad>
template <typename Scalar>
void
GenericClusterDynamicsGroupedNodalKernelTempl<is_ad>::fillConcentrationCache(
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
GenericClusterDynamicsGroupedNodalKernelTempl<is_ad>::fillFluxCache(
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
GenericClusterDynamicsGroupedNodalKernelTempl<is_ad>::computeGroupedResidualFromState(
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
GenericClusterDynamicsGroupedNodalKernelTempl<is_ad>::jacobianSetup()
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
GenericClusterDynamicsGroupedNodalKernelTempl<is_ad>::computeQpResidual(
    GenericRealEigenVector<is_ad> & residual)
{
  const auto n_comp = _u[_qp].size();
  if (n_comp != _layout.componentCount())
    mooseError("ClusterDynamicsGroupedNodalKernel expected ",
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
GenericClusterDynamicsGroupedNodalKernelTempl<false>::computeQpJacobian()
{
  const auto n_comp = static_cast<unsigned int>(_u[_qp].size());
  if (n_comp != _layout.componentCount())
    mooseError("ClusterDynamicsGroupedNodalKernel expected ",
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
GenericClusterDynamicsGroupedNodalKernelTempl<true>::computeQpJacobian()
{
  mooseError("Internal error: computeQpJacobian should never be called for the AD version");
}

// Non-AD off-diagonal Jacobian against the coupled vacancy array's own monomer component
// (component 0). Unlike computeQpJacobian()'s intra-variable derivative, concentration[] here
// does NOT depend on the coupled vacancy variable (only beta/alpha do, through
// monomerDiffusivity()), so there is no a0/a1 reconstruction chain rule to carry: every flux's
// derivative is simply fluxPrime(n) = betaPrime(n)*c1*concentration[n] -
// alphaPrime(n+1)*concentration[n+1], and every residual row is the exact same linear
// combination of fluxPrime(n) that computeQpJacobian() already uses for flux(n) -- written to
// the single column 0 of the coupled variable instead of fanned out across this kernel's own
// three affected columns.
template <>
void
GenericClusterDynamicsGroupedNodalKernelTempl<false>::computeQpOffDiagJacobian(
    const unsigned int jvar)
{
  if (!_vacancy_coupling_active || jvar != _vacancy_var_num)
    return;

  const auto n_comp = static_cast<unsigned int>(_u[_qp].size());
  ensureCoefficientCache();

  std::vector<Real> state(n_comp);
  for (unsigned int i = 0; i < n_comp; ++i)
    state[i] = _u[_qp](i);
  std::vector<Real> concentration;
  fillConcentrationCache(state, concentration);
  const Real c1 = concentration[1];

  std::vector<Real> flux_prime(_num_cluster_sizes + 1, 0.0);
  for (unsigned int n = 1; n < _num_cluster_sizes; ++n)
    flux_prime[n] = betaPrime(n) * c1 * concentration[n] - alphaPrime(n + 1) * concentration[n + 1];

  Real total = 0.0;
  auto addFluxPrimeDerivative = [&](unsigned int n, Real weight)
  {
    // J_0 and J_N vanish at the ends of the truncated size space, mirroring addFluxDerivative()
    if (n == 0 || n >= _num_cluster_sizes)
      return;
    total += weight * flux_prime[n];
  };

  // Monomer: same weights as computeQpJacobian()'s own monomer row
  total = 0.0;
  for (unsigned int n = 1; n < _num_cluster_sizes; ++n)
    addFluxPrimeDerivative(n, n == 1 ? 2.0 : 1.0);
  setJacobian(_layout.monomerComponent(), 0, total);

  // Explicit sizes
  for (unsigned int n = 2; n <= _layout.explicitMax(); ++n)
  {
    total = 0.0;
    addFluxPrimeDerivative(n - 1, -1.0);
    addFluxPrimeDerivative(n, 1.0);
    setJacobian(_layout.explicitComponent(n), 0, total);
  }

  for (const auto & bin : _layout.groups())
  {
    total = 0.0;
    addFluxPrimeDerivative(bin.start - 1, -1.0 / bin.width);
    addFluxPrimeDerivative(bin.end, 1.0 / bin.width);
    setJacobian(bin.l0_component, 0, total);

    if (bin.sigma2 > 0.0 && bin.width > 1.0)
    {
      const Real factor = (bin.width - 1.0) / (2.0 * bin.sigma2 * bin.width);
      total = 0.0;
      addFluxPrimeDerivative(bin.start - 1, factor);
      for (unsigned int n = bin.start; n < bin.end; ++n)
        addFluxPrimeDerivative(n, -2.0 * factor / (bin.width - 1.0));
      addFluxPrimeDerivative(bin.end, factor);
      setJacobian(bin.l1_component, 0, total);
    }
  }
}

template <>
void
GenericClusterDynamicsGroupedNodalKernelTempl<true>::computeQpOffDiagJacobian(const unsigned int)
{
  mooseError(
      "Internal error: computeQpOffDiagJacobian should never be called for the AD version");
}

template class GenericClusterDynamicsGroupedNodalKernelTempl<false>;
template class GenericClusterDynamicsGroupedNodalKernelTempl<true>;
