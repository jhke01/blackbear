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

#include "ContinuumClusterDynamicsNodalKernel.h"

#include "ClusterDynamicsSolverCheck.h"
#include "FEProblemBase.h"

#include "libmesh/utility.h"

registerMooseObject("BlackBearApp", ContinuumClusterDynamicsNodalKernel);

namespace
{
// The constants of GenericClusterDynamicsGroupedNodalKernel, so the coefficients agree with it
constexpr Real solute_kB = 1.380649e-23;
constexpr Real solute_eV = 1.602176634e-19;
} // namespace

InputParameters
ContinuumClusterDynamicsNodalKernel::validParams()
{
  InputParameters params = ArrayNodalKernel::validParams();
  params.addClassDescription(
      "Solute clustering with a size distribution that is discrete up to switch_size and a "
      "Fokker-Planck finite-volume continuum above it, the alternative to "
      "ClusterDynamicsGroupedNodalKernel (rate_model simple or interfacial_energy).");

  params.addParam<Real>("generation", 0.0, "Monomer generation rate G_1.");
  params.addParam<Real>("sink", 0.0, "Linear sink coefficient k_s of the monomer.");
  MooseEnum rate_model("simple interfacial_energy", "simple");
  params.addParam<MooseEnum>("rate_model",
                             rate_model,
                             "'simple': beta_n = beta0 n^(1/3), alpha_n = alpha0 n^(1/3). "
                             "'interfacial_energy': coefficients from the cluster geometry, the "
                             "interfacial energy and detailed balance.");
  params.addParam<Real>("beta0", 0.0, "Base absorption coefficient for rate_model = simple.");
  params.addParam<Real>("alpha0", 0.0, "Base emission coefficient for rate_model = simple.");
  params.addParam<Real>("temperature", 0.0, "Temperature T [K] for the interfacial-energy model.");
  MooseEnum diffusivity_model("constant arrhenius", "constant");
  params.addParam<MooseEnum>("diffusivity_model", diffusivity_model, "constant or arrhenius.");
  params.addParam<Real>("monomer_diffusivity", 0.0, "Monomer diffusivity [m^2/s] (constant).");
  params.addParam<Real>("D0", 0.0, "Diffusion prefactor [m^2/s] (arrhenius).");
  params.addParam<Real>("Q_eV", 0.0, "Activation energy [eV] (arrhenius).");
  params.addParam<Real>("radiation_enhanced_factor",
                        1.0,
                        "Multiplicative factor applied to the monomer diffusivity.");
  params.addCoupledVar("vacancy_coupled_var",
                       0.0,
                       "Optional vacancy array variable: the monomer diffusivity is multiplied by "
                       "its monomer concentration over the thermal value. It must have the "
                       "monomer in component 0.");
  params.addParam<Real>("vacancy_formation_energy_eV",
                        0.0,
                        "Vacancy formation energy [eV] for the thermal concentration "
                        "4 exp(-E_f / kT), required with vacancy_coupled_var.");
  MooseEnum binding_energy_model("interfacial_energy capillary", "interfacial_energy");
  params.addParam<MooseEnum>("binding_energy_model",
                             binding_energy_model,
                             "interfacial_energy: E_b = Omega - T DeltaS - 36 pi^(1/3) ... "
                             "(see ClusterDynamicsGroupedNodalKernel); capillary: A - B (n^p - "
                             "(n - 1)^p).");
  params.addParam<Real>("sigma", 0.0, "Interfacial energy [J/m^2].");
  params.addParam<Real>("atomic_volume", 0.0, "Atomic volume [m^3].");
  params.addParam<Real>("Omega_kB_K", 0.0, "Enthalpy term Omega [k_B K].");
  params.addParam<Real>("DeltaS_kB", 0.0, "Entropy term DeltaS [k_B].");
  params.addParam<Real>("binding_energy_constant_eV", 0.0, "Constant A [eV] (capillary).");
  params.addParam<Real>("binding_energy_enthalpy_kB_K", 0.0, "Enthalpy term H [k_B K] (capillary).");
  params.addParam<Real>("binding_energy_entropy_kB", 0.0, "Entropy term S [k_B] (capillary).");
  params.addParam<Real>("binding_energy_coefficient_eV", 0.0, "Coefficient B [eV] (capillary).");
  params.addParam<Real>("binding_energy_exponent", 2.0 / 3.0, "Exponent p (capillary).");

  params.addRequiredParam<unsigned int>("num_cluster_sizes", "Largest cluster size, N.");
  params.addRequiredParam<unsigned int>(
      "switch_size", "Largest cluster size that is an individual unknown; larger sizes are cells.");
  params.addParam<Real>("continuum_initial_width", 1.0, "Width of the first continuum cell.");
  params.addParam<Real>("continuum_width_ratio", 1.04, "Growth factor of the cell width.");
  params.addParam<bool>("use_limiter", true, "Minmod-limited second-order face values.");
  return params;
}

ContinuumClusterDynamicsNodalKernel::ContinuumClusterDynamicsNodalKernel(
    const InputParameters & parameters)
  : ArrayNodalKernel(parameters),
    _generation(getParam<Real>("generation")),
    _sink(getParam<Real>("sink")),
    _temperature(getParam<Real>("temperature")),
    _atomic_volume(getParam<Real>("atomic_volume")),
    _static_diffusivity(
        getParam<Real>("radiation_enhanced_factor") *
        (getParam<MooseEnum>("diffusivity_model") == "constant"
             ? getParam<Real>("monomer_diffusivity")
             : getParam<Real>("D0") *
                   std::exp(-getParam<Real>("Q_eV") * solute_eV / (solute_kB * _temperature)))),
    _simple(getParam<MooseEnum>("rate_model") == "simple"),
    _beta0(getParam<Real>("beta0")),
    _alpha0(getParam<Real>("alpha0")),
    _sigma(getParam<Real>("sigma")),
    _omega(getParam<Real>("Omega_kB_K") * solute_kB),
    _delta_s(getParam<Real>("DeltaS_kB") * solute_kB),
    _binding_constant(getParam<Real>("binding_energy_constant_eV") * solute_eV),
    _binding_enthalpy(getParam<Real>("binding_energy_enthalpy_kB_K") * solute_kB),
    _binding_entropy(getParam<Real>("binding_energy_entropy_kB") * solute_kB),
    _binding_coefficient(getParam<Real>("binding_energy_coefficient_eV") * solute_eV),
    _binding_exponent(getParam<Real>("binding_energy_exponent")),
    _capillary(getParam<MooseEnum>("binding_energy_model") == "capillary"),
    _vacancy_array(coupledArrayValue("vacancy_coupled_var")),
    _vacancy_active(isCoupled("vacancy_coupled_var")),
    _vacancy_var_num(coupled("vacancy_coupled_var")),
    _vacancy_c_thermal(
        _temperature > 0.0
            ? 4.0 * std::exp(-getParam<Real>("vacancy_formation_energy_eV") * solute_eV /
                             (solute_kB * _temperature))
            : 0.0),
    _layout(getParam<unsigned int>("num_cluster_sizes"),
            getParam<unsigned int>("switch_size"),
            getParam<Real>("continuum_initial_width"),
            getParam<Real>("continuum_width_ratio")),
    _flux(_layout, getParam<bool>("use_limiter"))
{
  if (!_simple)
  {
    if (_temperature <= 0.0)
      paramError("temperature", "must be positive for rate_model = interfacial_energy.");
    if (_atomic_volume <= 0.0)
      paramError("atomic_volume", "must be positive for rate_model = interfacial_energy.");
    if (_static_diffusivity <= 0.0)
      mooseError(type(), " requires a positive monomer diffusivity for rate_model = interfacial_energy.");
  }
  if (_vacancy_active)
  {
    if (_vacancy_var_num == _var.number())
      paramError("vacancy_coupled_var", "must be a different variable.");
    if (_simple)
      paramError("vacancy_coupled_var",
                 "has no effect with rate_model = simple, which does not use the diffusivity.");
    if (_vacancy_c_thermal <= 0.0)
      paramError("vacancy_formation_energy_eV", "is required with vacancy_coupled_var.");
  }
  if (_count != _layout.componentCount())
    mooseError(type(),
               " expects ",
               _layout.componentCount(),
               " array components from the layout, but the variable has ",
               _count,
               ".");
  if (!_fe_problem.useHashTableMatrixAssembly())
    mooseError(type(), " requires Problem/use_hash_table_matrix_assembly = true");
  if (_tid == 0)
    ClusterDynamics::checkSolverSetup(
        *this, _fe_problem, _var.sys().number(), _var.number(), _var.count());
}

Real
ContinuumClusterDynamicsNodalKernel::monomerDiffusivity() const
{
  Real diffusivity = _static_diffusivity;
  if (_vacancy_active)
    diffusivity *= _vacancy_array[_qp](0) / _vacancy_c_thermal;
  return diffusivity;
}

Real
ContinuumClusterDynamicsNodalKernel::radius(unsigned int n) const
{
  return std::cbrt(3.0 * static_cast<Real>(n) * _atomic_volume / (4.0 * libMesh::pi));
}

Real
ContinuumClusterDynamicsNodalKernel::bindingEnergy(unsigned int n) const
{
  const Real nn = static_cast<Real>(n);
  if (_capillary)
  {
    const Real constant_term =
        _binding_constant != 0.0 ? _binding_constant : _binding_enthalpy - _temperature * _binding_entropy;
    const Real coefficient = _binding_coefficient != 0.0
                                 ? _binding_coefficient
                                 : std::cbrt(36.0 * libMesh::pi) * std::pow(_atomic_volume, 2.0 / 3.0) * _sigma;
    return constant_term - coefficient * (std::pow(nn, _binding_exponent) - std::pow(nn - 1.0, _binding_exponent));
  }
  const Real surface = std::cbrt(36.0 * libMesh::pi) * std::pow(_atomic_volume, 2.0 / 3.0) * _sigma;
  return _omega - _temperature * _delta_s -
         surface * (std::pow(nn, 2.0 / 3.0) - std::pow(nn - 1.0, 2.0 / 3.0));
}

Real
ContinuumClusterDynamicsNodalKernel::betaShape(unsigned int n) const
{
  if (n == 0)
    return 0.0;
  if (_simple)
    return _beta0 * std::cbrt(static_cast<Real>(n));
  return 4.0 * libMesh::pi * (radius(1) + radius(n)) / _atomic_volume;
}

Real
ContinuumClusterDynamicsNodalKernel::alphaShape(unsigned int n) const
{
  if (n < 2)
    return 0.0;
  if (_simple)
    return _alpha0 * std::cbrt(static_cast<Real>(n));
  // alpha(n) = beta(n - 1) exp(-E_b(n) / kT), the same detailed balance as the grouped kernel
  return betaShape(n - 1) * std::exp(-bindingEnergy(n) / (solute_kB * _temperature));
}

void
ContinuumClusterDynamicsNodalKernel::ensureShapes()
{
  if (_have_shapes)
    return;
  _shapes.build(
      _layout,
      [this](unsigned int n) { return betaShape(n); },
      [this](unsigned int n) { return alphaShape(n); });
  _coefficients = _shapes;
  _have_shapes = true;
  _applied_diffusivity = _simple ? 1.0 : -1.0;
}

void
ContinuumClusterDynamicsNodalKernel::applyDiffusivity()
{
  ensureShapes();
  if (_simple)
    return;
  // beta and alpha are proportional to the diffusivity, so the tables are the unit-diffusivity
  // shapes times the current value, rebuilt only when it changed
  const Real diffusivity = monomerDiffusivity();
  if (diffusivity == _applied_diffusivity)
    return;
  auto scale = [diffusivity](const std::vector<Real> & from, std::vector<Real> & to)
  {
    to.resize(from.size());
    for (const auto i : index_range(from))
      to[i] = from[i] * diffusivity;
  };
  scale(_shapes.beta_discrete, _coefficients.beta_discrete);
  scale(_shapes.alpha_discrete, _coefficients.alpha_discrete);
  scale(_shapes.beta_center, _coefficients.beta_center);
  scale(_shapes.alpha_center, _coefficients.alpha_center);
  scale(_shapes.beta_face, _coefficients.beta_face);
  scale(_shapes.alpha_face, _coefficients.alpha_face);
  _applied_diffusivity = diffusivity;
}

void
ContinuumClusterDynamicsNodalKernel::computeQpResidual(RealEigenVector & residual)
{
  applyDiffusivity();
  _flux.rate(&_u[_qp](0), _coefficients, _generation, _sink, _rate);
  residual.resize(_count);
  for (const auto i : make_range(_count))
    residual(i) = -_rate[i];
}

void
ContinuumClusterDynamicsNodalKernel::computeQpJacobian()
{
  applyDiffusivity();
  _flux.jacobian(&_u[_qp](0), _coefficients, _sink, _triplets);
  for (const auto & t : _triplets)
    setJacobian(t.row, t.column, t.value);
}

void
ContinuumClusterDynamicsNodalKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (!_vacancy_active || jvar != _vacancy_var_num)
    return;

  // The rates are linear in the diffusivity, which is linear in the vacancy monomer C_v:
  // d(rate) / d(C_v) = rate_diffusivity_part / C_v. The sink and the generation do not depend on
  // it, so the derivative is that of the rate evaluated at unit generation and sink zero.
  applyDiffusivity();
  std::vector<Real> rate_without_sources;
  _flux.rate(&_u[_qp](0), _coefficients, 0.0, 0.0, rate_without_sources);
  const Real d_scale = _static_diffusivity / _vacancy_c_thermal / monomerDiffusivity();
  for (const auto i : make_range(_count))
    setJacobian(i, 0, -rate_without_sources[i] * d_scale);
}
