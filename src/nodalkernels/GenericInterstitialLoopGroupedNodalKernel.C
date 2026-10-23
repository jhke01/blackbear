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
  InputParameters params = GenericDefectFamilyGroupedNodalKernel<
      GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>, is_ad>::validParams();
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
  params.addParam<Real>(
      "z_i", 1.0, "Large-loop limit Z_i of the dislocation bias, approached as n -> infinity.");
  params.addParam<Real>(
      "z_loop_i",
      1.0,
      "Small-loop coefficient Z_loop_i: the bias at n = 1 is sqrt(b/(8 pi a0)) * Z_loop_i.");
  params.addParam<Real>(
      "alpha_loop_i", 1.0, "Exponent controlling how fast Z_i(n) approaches Z_i with size.");
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

  MooseEnum binding_energy_model("capillary binding_energy_table formation_energy_table",
                                 "capillary");
  params.addParam<MooseEnum>(
      "binding_energy_model",
      binding_energy_model,
      "Binding-energy model for loop emission. 'capillary' (default) always uses the "
      "loop_energy_coefficient_eV/loop_energy_constant_eV law above. 'binding_energy_table' or "
      "'formation_energy_table' use direct DFT/MD data for small n, falling back to the "
      "capillary law for any n beyond the table's range -- the table's own length is the "
      "crossover point.");
  params.addParam<std::vector<Real>>(
      "binding_energy_table_eV",
      std::vector<Real>(),
      "Direct binding energies [eV] for binding_energy_model = binding_energy_table. Entry 0 "
      "corresponds to cluster size n = 2.");
  params.addParam<std::vector<Real>>(
      "formation_energy_table_eV",
      std::vector<Real>(),
      "Formation energies [eV] for binding_energy_model = formation_energy_table. Entry 0 "
      "corresponds to cluster size n = 1.");

  params.addParam<Real>(
      "D0_v", 0.0, "Vacancy diffusion prefactor D0_v [m^2/s] (opposite species), for this "
      "loop's own absorption of a vacancy monomer.");
  params.addParam<Real>(
      "Em_v_eV", 0.0, "Vacancy migration energy Em_v [eV] (opposite species).");
  params.addParam<Real>(
      "z_v",
      1.0,
      "Large-loop limit Z_v of the bias for loops absorbing vacancies (n -> infinity).");
  params.addParam<Real>(
      "z_loop_v",
      1.0,
      "Small-loop coefficient Z_loop_v for loops absorbing vacancies: the bias at n = 1 is "
      "sqrt(b/(8 pi a0)) * Z_loop_v.");
  params.addParam<Real>(
      "alpha_loop_v",
      1.0,
      "Exponent controlling how fast the loop vacancy-absorption bias approaches Z_v with "
      "size.");

  MooseEnum opposite_geometry("void loop", "void");
  params.addParam<MooseEnum>(
      "opposite_geometry",
      opposite_geometry,
      "Geometry of the coupled opposite family (opposite_coupled_var): 'void' (default, a "
      "three-dimensional cluster with the finite-sink-strength absorption coefficient) or 'loop' "
      "(a two-dimensional cluster, e.g. a vacancy modeled as a loop following Pokor et al. 2004). "
      "This selects how the opposite family's absorption of THIS kernel's monomer is summed for "
      "the monomer-pool loss (Piece B), and must match how the opposite kernel computes its own "
      "absorption (Piece A).");
  params.addParam<Real>("opposite_burgers_vector",
                        0.0,
                        "Burgers vector b [m] of the opposite family's loops "
                        "(opposite_geometry = loop only). Defaults to this kernel's own "
                        "burgers_vector when left at zero.");
  params.addParam<Real>("opposite_z_i",
                        1.0,
                        "Baseline bias Z of the opposite family's small loops for capturing this "
                        "kernel's monomer (opposite_geometry = loop only; required to be set).");
  params.addParam<Real>("opposite_z_loop_i",
                        1.0,
                        "Small-loop coefficient Z_loop of the opposite family's loops for "
                        "capturing this kernel's monomer: their bias at n = 1 is "
                        "sqrt(b/(8 pi a0)) * Z_loop (opposite_geometry = loop only; required to be "
                        "set).");
  params.addParam<Real>("opposite_alpha_loop_i",
                        1.0,
                        "Exponent controlling how fast the opposite family's loop bias approaches "
                        "opposite_z_i with size (opposite_geometry = loop only; required to be "
                        "set).");
  return params;
}

template <bool is_ad>
GenericInterstitialLoopGroupedNodalKernelTempl<
    is_ad>::GenericInterstitialLoopGroupedNodalKernelTempl(const InputParameters & parameters)
  : GenericDefectFamilyGroupedNodalKernel<GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>,
                                          is_ad>(parameters),
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
    _binding_energy_model(this->template getParam<MooseEnum>("binding_energy_model")
                              .template getEnum<BindingEnergyModel>()),
    _binding_energy_table_eV(this->template getParam<std::vector<Real>>("binding_energy_table_eV")),
    _formation_energy_table_eV(
        this->template getParam<std::vector<Real>>("formation_energy_table_eV")),
    _opposite_D0(this->template getParam<Real>("D0_v")),
    _opposite_Q(this->template getParam<Real>("Em_v_eV") * loop_eV),
    _z_v(this->template getParam<Real>("z_v")),
    _z_loop_v(this->template getParam<Real>("z_loop_v")),
    _alpha_loop_v(this->template getParam<Real>("alpha_loop_v")),
    _opposite_is_loop(this->template getParam<MooseEnum>("opposite_geometry") == "loop"),
    _opposite_burgers_vector(this->template getParam<Real>("opposite_burgers_vector") != 0.0
                                 ? this->template getParam<Real>("opposite_burgers_vector")
                                 : _burgers_vector),
    _opposite_z_i(this->template getParam<Real>("opposite_z_i")),
    _opposite_z_loop_i(this->template getParam<Real>("opposite_z_loop_i")),
    _opposite_alpha_loop_i(this->template getParam<Real>("opposite_alpha_loop_i"))
{
  if (_temperature <= 0.0)
    mooseError(this->type(), " requires temperature > 0.");
  if (_D0 <= 0.0)
    mooseError(this->type(), " requires D0_i > 0.");
  if (_Q <= 0.0)
    mooseError(this->type(), " requires Em_i_eV > 0.");
  if (_opposite_burgers_vector <= 0.0)
    this->paramError("opposite_burgers_vector", "must be positive, or zero to use burgers_vector.");
  if (_atomic_volume <= 0.0)
    mooseError(this->type(), " requires atomic_volume > 0.");
  if (_burgers_vector <= 0.0)
    mooseError(this->type(), " requires lattice_parameter > 0 or burgers_vector > 0.");
  if (_lattice_parameter <= 0.0)
    mooseError(this->type(), " requires lattice_parameter > 0.");
  if (_formation_energy_eV <= 0.0)
    mooseError(this->type(), " requires formation_energy_eV > 0.");
  if (_energy_exponent <= 0.0)
    mooseError(this->type(), " requires loop_energy_exponent > 0.");
  if (_energy_coefficient_eV <= 0.0)
    mooseError(this->type(), " requires loop_energy_coefficient_eV > 0.");
  if (_binding_energy_model == BindingEnergyModel::BINDING_ENERGY_TABLE &&
      _binding_energy_table_eV.empty())
    mooseError(this->type(), " with binding_energy_model = binding_energy_table requires "
                             "binding_energy_table_eV.");
  if (_binding_energy_model == BindingEnergyModel::FORMATION_ENERGY_TABLE &&
      _formation_energy_table_eV.size() < 2)
    mooseError(this->type(), " with binding_energy_model = formation_energy_table requires at "
                             "least E_f(1) and E_f(2).");
  // The opposite-absorption parameters only need to be physically valid when the opposite
  // coupling is actually active (see the base class constructor for why).
  if (_opposite_coupling_active)
  {
    if (_opposite_D0 <= 0.0)
      mooseError(this->type(), " requires D0_v > 0 when opposite_coupled_var is provided.");
    if (_opposite_Q <= 0.0)
      mooseError(this->type(), " requires Em_v_eV > 0 when opposite_coupled_var is provided.");
  }
  if (_opposite_is_loop && _opposite_coupling_active)
  {
    // A loop partner's bias is the partner kernel's own Piece A bias for capturing this kernel's
    // monomer, so there is no sensible default: silently using 1.0 would give a wrong coefficient.
    for (const auto & name : {"opposite_z_i", "opposite_z_loop_i", "opposite_alpha_loop_i"})
      if (!this->parameters().isParamSetByUser(name))
        this->paramError(name,
                         "must be set when opposite_geometry = loop: it is the opposite loop "
                         "family's bias for capturing this kernel's monomer.");
  }
  if (!_opposite_is_loop && (this->parameters().isParamSetByUser("opposite_z_i") ||
                             this->parameters().isParamSetByUser("opposite_z_loop_i") ||
                             this->parameters().isParamSetByUser("opposite_alpha_loop_i") ||
                             this->parameters().isParamSetByUser("opposite_burgers_vector")))
    this->paramWarning("opposite_geometry",
                       "opposite_z_i, opposite_z_loop_i, opposite_alpha_loop_i and "
                       "opposite_burgers_vector have no effect unless opposite_geometry = loop.");
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
  // Z_i(n) interpolates from sqrt(b/(8 pi a0)) * z_loop_i at n = 1 to the large-loop bias z_i as n
  // grows,
  // following the dislocation-loop capture-efficiency model of Wiedersich/Jourdan-style AMSS
  // rate theory (DD_AMSS_Z_I, DD_AMSS_Z_LOOP_I, DD_AMSS_ALPHA_LOOP_I in the reference model).
  Real z = std::sqrt(_burgers_vector / (8.0 * libMesh::pi * _lattice_parameter)) * _z_loop_i - _z_i;
  z *= std::pow(1.0 / static_cast<Real>(n), 0.5 * _alpha_loop_i);
  return z + _z_i;
}

template <bool is_ad>
Real
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::bindingEnergy(unsigned int n) const
{
  if (n <= 1)
    mooseError(this->type(), " requires n >= 2 for bindingEnergy(n).");

  if (_binding_energy_model == BindingEnergyModel::BINDING_ENERGY_TABLE &&
      n - 2 < _binding_energy_table_eV.size())
    return _binding_energy_table_eV[n - 2];

  if (_binding_energy_model == BindingEnergyModel::FORMATION_ENERGY_TABLE &&
      n <= _formation_energy_table_eV.size())
    return _formation_energy_table_eV[n - 2] + _formation_energy_table_eV[0] -
           _formation_energy_table_eV[n - 1];

  // Capillary law: the default mode, and the fallback beyond either table's range.
  return _formation_energy_eV - (formationEnergy(n) - formationEnergy(n - 1));
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

  return beta(n - 1) * std::exp(-bindingEnergy(n) * loop_eV / (loop_kB * _temperature));
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
  if (_opposite_is_loop)
    return std::sqrt(static_cast<Real>(n) * _atomic_volume /
                     (libMesh::pi * _opposite_burgers_vector));
  return std::cbrt(3.0 * static_cast<Real>(n) * _atomic_volume / (4.0 * libMesh::pi));
}

template <bool is_ad>
Real
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::oppositeLoopBias(unsigned int n) const
{
  // Same functional form as loopBias(), with the partner's Burgers vector and bias parameters: this
  // is what the partner kernel evaluates as its own biasVacancy() when it absorbs this kernel's
  // monomer, so the two sides of one absorption event use one coefficient.
  Real z = std::sqrt(_opposite_burgers_vector / (8.0 * libMesh::pi * _lattice_parameter)) *
               _opposite_z_loop_i -
           _opposite_z_i;
  z *= std::pow(1.0 / static_cast<Real>(n), 0.5 * _opposite_alpha_loop_i);
  return z + _opposite_z_i;
}

template <bool is_ad>
Real
GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>::oppositeBetaOpposite(unsigned int n) const
{
  if (n < _opposite_beta_opposite_cache.size())
    return _opposite_beta_opposite_cache[n];

  const Real r = oppositeRadius(n);
  if (_opposite_is_loop)
    return 2.0 * libMesh::pi * r * oppositeLoopBias(n) * monomerDiffusivity() / _atomic_volume;
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
        _opposite_is_loop
            ? 2.0 * libMesh::pi * r * oppositeLoopBias(n) * diffusivity / _atomic_volume
            : 4.0 * libMesh::pi * r * diffusivity /
                  ((1.0 + _lattice_parameter / r) * _atomic_volume);
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
    alpha_values[n] =
        beta_values[n - 1] * std::exp(-bindingEnergy(n) * loop_eV / (loop_kB * _temperature));

  _beta_cache = std::move(beta_values);
  _alpha_cache = std::move(alpha_values);
  _beta_opposite_cache = std::move(beta_opposite_values);
  _cache_size = _num_cluster_sizes;
}

template class GenericInterstitialLoopGroupedNodalKernelTempl<false>;
template class GenericInterstitialLoopGroupedNodalKernelTempl<true>;
