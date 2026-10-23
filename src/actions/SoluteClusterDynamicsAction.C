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

#include "SoluteClusterDynamicsAction.h"

#include <cmath>
#include "DefectSpeciesAction.h"
#include "DefectDynamicsAction.h"

#include "FEProblemBase.h"

registerMooseAction("BlackBearApp", SoluteClusterDynamicsAction, "add_variable");
registerMooseAction("BlackBearApp", SoluteClusterDynamicsAction, "add_function");
registerMooseAction("BlackBearApp", SoluteClusterDynamicsAction, "add_ic");
registerMooseAction("BlackBearApp", SoluteClusterDynamicsAction, "add_nodal_kernel");
registerMooseAction("BlackBearApp", SoluteClusterDynamicsAction, "add_postprocessor");

InputParameters
SoluteClusterDynamicsAction::validParams()
{
  InputParameters params = Action::validParams();
  MooseEnum rate_model("simple interfacial_energy", "simple");
  MooseEnum diffusivity_model("constant arrhenius", "constant");
  MooseEnum binding_energy_model("interfacial_energy capillary binding_energy_table "
                                 "formation_energy_table",
                                 "interfacial_energy");
  MooseEnum group_spacing("linear geometric", "linear");
  params.addClassDescription(
      "A front end to ClusterDynamicsGroupedNodalKernel, the solute-clustering kernel (used for "
      "Cu precipitation, for example): builds the solute array variable, initial condition, "
      "time-derivative/dynamics kernel, and monomer/total-density postprocessors from the "
      "physical parameters below instead of requiring [GroupedVariables]/[NodalKernels] to be "
      "hand-typed. When a [DefectDynamics] block is also present in the same input, this action "
      "automatically couples the solute monomer diffusivity to that block's vacancy species "
      "(see vacancy_coupled_var on GenericClusterDynamicsGroupedNodalKernel) -- no separate "
      "switch is needed; declaring both blocks together activates concurrent point-defect "
      "clustering and its effect on solute diffusion.");

  params.addRequiredParam<Real>(
      "initial_concentration", "Initial solute monomer concentration (atomic fraction).");

  params.addParam<Real>("generation", 0.0, "Monomer generation rate G_1 (leave at zero for "
                        "solute clustering such as Cu precipitation).");
  params.addParam<Real>("sink", 0.0, "Linear sink coefficient k_s for monomer loss (leave at "
                        "zero for solute clustering such as Cu precipitation).");

  params.addParam<MooseEnum>(
      "tail",
      MooseEnum("grouping fokker_planck", "grouping"),
      "How the large cluster sizes are stored. 'grouping' (default): bins with a linear "
      "reconstruction, set by explicit_max and the group_* parameters. 'fokker_planck': sizes up to "
      "switch_size are individual unknowns and larger sizes are finite-volume cells of a "
      "Fokker-Planck continuum (Pokor et al. 2004, Section 4.1), set by switch_size and the "
      "continuum_* parameters. Supports the interfacial_energy and capillary binding-energy "
      "models, not the tables.");
  params.addParam<unsigned int>("switch_size",
                                200,
                                "Largest cluster size that is an individual unknown when tail = "
                                "fokker_planck; larger sizes are continuum cells.");
  params.addParam<Real>("continuum_initial_width",
                        1.0,
                        "Width of the first continuum cell when tail = fokker_planck.");
  params.addParam<Real>("continuum_width_ratio",
                        1.04,
                        "Growth factor of the cell width from one cell to the next when tail = "
                        "fokker_planck.");
  params.addParam<bool>("continuum_limiter",
                        true,
                        "Reconstruct the cell values at a face with a minmod-limited slope "
                        "(second order) when tail = fokker_planck; if false the cell averages are "
                        "used (first-order upwind).");

  params.addRequiredParam<unsigned int>(
      "num_cluster_sizes", "Maximum physical cluster size represented by the grouped array.");
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
  params.addParam<Real>("beta0", 0.0, "Base absorption coefficient for rate_model = simple.");
  params.addParam<Real>("alpha0", 0.0, "Base emission coefficient for rate_model = simple.");

  params.addParam<MooseEnum>(
      "diffusivity_model",
      diffusivity_model,
      "Diffusivity model used only when rate_model = interfacial_energy. 'constant' uses "
      "monomer_diffusivity directly, while 'arrhenius' computes the diffusivity from D0, Q_eV, "
      "and temperature.");
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

  params.addParam<MooseEnum>(
      "binding_energy_model",
      binding_energy_model,
      "Binding-energy model used by detailed-balance emission when rate_model = "
      "interfacial_energy. 'interfacial_energy' uses Omega, DeltaS, sigma, and atomic volume; "
      "'capillary' uses A - B*(n^p - (n-1)^p); 'binding_energy_table' reads direct binding "
      "energies for n >= 2; and 'formation_energy_table' computes "
      "E_b(n) = E_f(n-1) + E_f(1) - E_f(n).");
  params.addParam<Real>(
      "sigma",
      0.0,
      "Interfacial energy sigma [J/m^2] for binding_energy_model = interfacial_energy or "
      "capillary when the capillary surface coefficient is computed from sigma.");
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
      "binding_energy_model = capillary.");
  params.addParam<Real>("binding_energy_entropy_kB",
                        0.0,
                        "Entropy-like term S [k_B] for the temperature-dependent constant in "
                        "binding_energy_model = capillary.");
  params.addParam<Real>(
      "binding_energy_coefficient_eV",
      0.0,
      "Direct coefficient B [eV] in binding_energy_model = capillary. Leave this at zero and "
      "provide sigma to compute B from interfacial energy instead.");
  params.addParam<Real>(
      "binding_energy_exponent", 2.0 / 3.0, "Exponent p in binding_energy_model = capillary.");
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
      "resolvable_radius",
      0.0,
      "If positive, adds a number-density and average-radius postprocessor pair "
      "(<name>_resolvable_density, <name>_resolvable_avg_radius) counting only clusters at or "
      "above this radius [m] -- the smallest cluster resolvable by whatever characterization "
      "technique the input is compared against (e.g. atom probe). Same name and meaning as "
      "resolvable_radius on a [DefectDynamics] species. The radius is converted to a cluster size "
      "with the same spherical law the solute kernel uses, r_n = (3 n V_at / 4 pi)^(1/3), "
      "using atomic_volume. Left at zero (the default), no such postprocessors are added.");
  params.addDeprecatedParam<unsigned int>(
      "density_n_minimum",
      1,
      "Smallest cluster size included in the <name>_total_density postprocessor.",
      "density_n_minimum is deprecated; <name>_total_density keeps counting every cluster size "
      "from this value up. Use resolvable_radius, which adds <name>_resolvable_density and "
      "<name>_resolvable_avg_radius.");

  return params;
}

SoluteClusterDynamicsAction::SoluteClusterDynamicsAction(const InputParameters & params)
  : Action(params)
{
}

bool
SoluteClusterDynamicsAction::usesContinuumTail() const
{
  return getParam<MooseEnum>("tail") == "fokker_planck";
}

ClusterContinuumLayout
SoluteClusterDynamicsAction::continuumLayout() const
{
  return ClusterContinuumLayout(getParam<unsigned int>("num_cluster_sizes"),
                                getParam<unsigned int>("switch_size"),
                                getParam<Real>("continuum_initial_width"),
                                getParam<Real>("continuum_width_ratio"));
}

void
SoluteClusterDynamicsAction::validateTail() const
{
  const char * continuum_params[] = {
      "switch_size", "continuum_initial_width", "continuum_width_ratio", "continuum_limiter"};
  if (!usesContinuumTail())
  {
    for (const auto * name : continuum_params)
      if (parameters().isParamSetByUser(name))
        paramWarning(name, name, " has no effect with tail = grouping.");
    return;
  }

  const char * grouping_params[] = {"use_grouping",
                                    "explicit_max",
                                    "group_spacing",
                                    "group_linear_width",
                                    "group_geometric_ratio",
                                    "enforce_group_nonnegative",
                                    "group_nonnegative_tolerance_factor"};
  for (const auto * name : grouping_params)
    if (parameters().isParamSetByUser(name))
      paramError(name, name, " only applies to tail = grouping.");

  const auto binding = getParam<MooseEnum>("binding_energy_model");
  if (getParam<MooseEnum>("rate_model") == "interfacial_energy" &&
      (binding == "binding_energy_table" || binding == "formation_energy_table"))
    paramError("binding_energy_model",
               "tail = fokker_planck supports only the interfacial_energy and capillary "
               "binding-energy models, not the tables.");
  if (getParam<unsigned int>("switch_size") < 2)
    paramError("switch_size", "must be at least 2.");
  if (getParam<unsigned int>("switch_size") >= getParam<unsigned int>("num_cluster_sizes"))
    paramError("switch_size", "must be smaller than num_cluster_sizes.");
  if (getParam<Real>("continuum_initial_width") <= 0.0)
    paramError("continuum_initial_width", "must be positive.");
  if (getParam<Real>("continuum_width_ratio") < 1.0)
    paramError("continuum_width_ratio", "must be at least 1.");

  // The resolvable population starts at the size n_min of the smallest resolvable radius. When n_min
  // lies above switch_size it falls inside a continuum cell and the resolvable density is only as
  // accurate as the cut of that cell.
  const Real resolvable_radius = getParam<Real>("resolvable_radius");
  const Real atomic_volume = getParam<Real>("atomic_volume");
  if (resolvable_radius > 0.0 && atomic_volume > 0.0)
  {
    const Real r1 = std::cbrt(3.0 * atomic_volume / (4.0 * libMesh::pi));
    const Real n_min = std::ceil(std::pow(resolvable_radius / r1, 3.0));
    if (static_cast<Real>(getParam<unsigned int>("switch_size")) < n_min)
      paramWarning("switch_size",
                   "is below the size ",
                   static_cast<unsigned int>(n_min),
                   " at which clusters become resolvable (resolvable_radius), so the cut of the "
                   "resolvable population falls inside a continuum cell and its density is less "
                   "accurate. Use switch_size of at least ",
                   static_cast<unsigned int>(n_min),
                   ".");
  }
}

std::string
SoluteClusterDynamicsAction::variableName() const
{
  return name() + "_clusters";
}

std::string
SoluteClusterDynamicsAction::functionName() const
{
  return name() + "_index_fn";
}

const DefectSpeciesAction *
SoluteClusterDynamicsAction::findVacancySpecies() const
{
  const auto defect_dynamics_blocks = _awh.getActions<DefectDynamicsAction>();
  if (defect_dynamics_blocks.empty())
    return nullptr;
  if (defect_dynamics_blocks.size() > 1)
    mooseError("[SoluteClusterDynamics] found ", defect_dynamics_blocks.size(), " [DefectDynamics] "
              "blocks; at most one is supported for automatic point-defect coupling.");

  const auto species = _awh.getActions<DefectSpeciesAction>();

  const DefectSpeciesAction * void_species = nullptr;
  const DefectSpeciesAction * named_vacancy = nullptr;
  for (const auto * sp : species)
  {
    if (sp->getParam<MooseEnum>("geometry") == "void")
    {
      if (void_species)
        mooseError("[SoluteClusterDynamics] found more than one geometry = void species in the "
                  "[DefectDynamics] block; the vacancy species to couple to is ambiguous.");
      void_species = sp;
    }
    if (sp->name() == "vacancy")
      named_vacancy = sp;
  }

  if (void_species)
    return void_species;
  // Fallback: a Pokor-style two-loop-species setup models vacancy as geometry = loop (see
  // DefectDynamicsAction's own two-loop-species support), so geometry alone cannot identify it
  // there; a species literally named "vacancy" is the fallback key.
  return named_vacancy;
}

void
SoluteClusterDynamicsAction::act()
{
  if (_current_task == "add_variable")
    validateTail();

  if (_current_task == "add_variable")
    actAddVariables();
  else if (_current_task == "add_function")
    actAddFunctions();
  else if (_current_task == "add_ic")
    actAddICs();
  else if (_current_task == "add_nodal_kernel")
    actAddNodalKernels();
  else if (_current_task == "add_postprocessor")
    actAddPostprocessors();
}

void
SoluteClusterDynamicsAction::actAddVariables()
{
  if (usesContinuumTail())
  {
    InputParameters params = _factory.getValidParams("ArrayMooseVariable");
    params.set<MooseEnum>("family") = "LAGRANGE";
    params.set<MooseEnum>("order") = "FIRST";
    params.set<unsigned int>("components") = continuumLayout().componentCount();
    _problem->addVariable("ArrayMooseVariable", variableName(), params);
    return;
  }

  const ClusterGroupingLayout layout(getParam<unsigned int>("num_cluster_sizes"),
                                     getParam<unsigned int>("explicit_max"),
                                     getParam<bool>("use_grouping"),
                                     getParam<MooseEnum>("group_spacing")
                                         .getEnum<ClusterGroupingSpacing>(),
                                     getParam<unsigned int>("group_linear_width"),
                                     getParam<Real>("group_geometric_ratio"));

  InputParameters params = _factory.getValidParams("ArrayMooseVariable");
  params.set<MooseEnum>("family") = "LAGRANGE";
  params.set<MooseEnum>("order") = "FIRST";
  params.set<unsigned int>("components") = layout.componentCount();
  _problem->addVariable("ArrayMooseVariable", variableName(), params);
}

void
SoluteClusterDynamicsAction::actAddFunctions()
{
  InputParameters params = _factory.getValidParams("PiecewiseConstant");
  params.set<std::vector<Real>>("x") = {0.0, 0.9};
  params.set<std::vector<Real>>("y") = {getParam<Real>("initial_concentration"), 0.0};
  params.set<MooseEnum>("direction") = "left";
  _problem->addFunction("PiecewiseConstant", functionName(), params);
}

void
SoluteClusterDynamicsAction::actAddICs()
{
  if (usesContinuumTail())
  {
    // Only the monomer is populated: component 0 reads the function at x = 0 and every other
    // component reads it at x >= 1, where it is zero
    InputParameters params = _factory.getValidParams("ArrayFunctionIndexIC");
    params.set<VariableName>("variable") = variableName();
    params.set<FunctionName>("function") = functionName();
    _problem->addInitialCondition("ArrayFunctionIndexIC", name() + "_ic", params);
    return;
  }

  InputParameters params = _factory.getValidParams("GroupedArrayFunctionIndexIC");
  params.set<VariableName>("variable") = variableName();
  params.set<FunctionName>("function") = functionName();
  params.set<unsigned int>("num_cluster_sizes") = getParam<unsigned int>("num_cluster_sizes");
  params.set<bool>("use_grouping") = getParam<bool>("use_grouping");
  params.set<unsigned int>("explicit_max") = getParam<unsigned int>("explicit_max");
  params.set<MooseEnum>("group_spacing") = getParam<MooseEnum>("group_spacing");
  params.set<unsigned int>("group_linear_width") = getParam<unsigned int>("group_linear_width");
  params.set<Real>("group_geometric_ratio") = getParam<Real>("group_geometric_ratio");
  _problem->addInitialCondition("GroupedArrayFunctionIndexIC", name() + "_ic", params);
}

void
SoluteClusterDynamicsAction::actAddNodalKernels()
{
  const std::string var_name = variableName();

  {
    InputParameters params = _factory.getValidParams("ArrayTimeDerivativeNodalKernel");
    params.set<NonlinearVariableName>("variable") = var_name;
    _problem->addNodalKernel("ArrayTimeDerivativeNodalKernel", name() + "_dot", params);
  }

  const std::string kernel_type =
      usesContinuumTail() ? "ContinuumClusterDynamicsNodalKernel" : "ClusterDynamicsGroupedNodalKernel";
  InputParameters params = _factory.getValidParams(kernel_type);
  params.set<NonlinearVariableName>("variable") = var_name;
  params.set<Real>("generation") = getParam<Real>("generation");
  params.set<Real>("sink") = getParam<Real>("sink");
  params.set<unsigned int>("num_cluster_sizes") = getParam<unsigned int>("num_cluster_sizes");
  if (usesContinuumTail())
  {
    params.set<unsigned int>("switch_size") = getParam<unsigned int>("switch_size");
    params.set<Real>("continuum_initial_width") = getParam<Real>("continuum_initial_width");
    params.set<Real>("continuum_width_ratio") = getParam<Real>("continuum_width_ratio");
    params.set<bool>("use_limiter") = getParam<bool>("continuum_limiter");
  }
  else
  {
    params.set<bool>("use_grouping") = getParam<bool>("use_grouping");
    params.set<unsigned int>("explicit_max") = getParam<unsigned int>("explicit_max");
    params.set<MooseEnum>("group_spacing") = getParam<MooseEnum>("group_spacing");
    params.set<unsigned int>("group_linear_width") = getParam<unsigned int>("group_linear_width");
    params.set<Real>("group_geometric_ratio") = getParam<Real>("group_geometric_ratio");
    params.set<bool>("enforce_group_nonnegative") = getParam<bool>("enforce_group_nonnegative");
    params.set<Real>("group_nonnegative_tolerance_factor") =
        getParam<Real>("group_nonnegative_tolerance_factor");
  }
  params.set<MooseEnum>("rate_model") = getParam<MooseEnum>("rate_model");
  params.set<Real>("beta0") = getParam<Real>("beta0");
  params.set<Real>("alpha0") = getParam<Real>("alpha0");
  params.set<MooseEnum>("diffusivity_model") = getParam<MooseEnum>("diffusivity_model");
  params.set<Real>("temperature") = getParam<Real>("temperature");
  params.set<Real>("monomer_diffusivity") = getParam<Real>("monomer_diffusivity");
  params.set<Real>("D0") = getParam<Real>("D0");
  params.set<Real>("Q_eV") = getParam<Real>("Q_eV");
  params.set<Real>("radiation_enhanced_factor") = getParam<Real>("radiation_enhanced_factor");
  // validateTail() has rejected the tables for the continuum, whose kernel has only the first two
  params.set<MooseEnum>("binding_energy_model") =
      std::string(getParam<MooseEnum>("binding_energy_model"));
  params.set<Real>("sigma") = getParam<Real>("sigma");
  params.set<Real>("atomic_volume") = getParam<Real>("atomic_volume");
  params.set<Real>("Omega_kB_K") = getParam<Real>("Omega_kB_K");
  params.set<Real>("DeltaS_kB") = getParam<Real>("DeltaS_kB");
  params.set<Real>("binding_energy_constant_eV") = getParam<Real>("binding_energy_constant_eV");
  params.set<Real>("binding_energy_enthalpy_kB_K") =
      getParam<Real>("binding_energy_enthalpy_kB_K");
  params.set<Real>("binding_energy_entropy_kB") = getParam<Real>("binding_energy_entropy_kB");
  params.set<Real>("binding_energy_coefficient_eV") =
      getParam<Real>("binding_energy_coefficient_eV");
  params.set<Real>("binding_energy_exponent") = getParam<Real>("binding_energy_exponent");
  if (!usesContinuumTail())
  {
    params.set<std::vector<Real>>("binding_energy_table_eV") =
        getParam<std::vector<Real>>("binding_energy_table_eV");
    params.set<std::vector<Real>>("formation_energy_table_eV") =
        getParam<std::vector<Real>>("formation_energy_table_eV");
  }

  const auto * vacancy = findVacancySpecies();
  if (vacancy)
  {
    // The vacancy supersaturation C_v / C_v_thermal divides the live vacancy concentration, which
    // [DefectDynamics] computes at its temperature, by the thermal value, which the solute kernel
    // evaluates at its own. Two temperatures would make the ratio wrong at equilibrium.
    const Real defect_temperature =
        _awh.getActions<DefectDynamicsAction>().front()->getParam<Real>("temperature");
    if (getParam<MooseEnum>("rate_model") == "interfacial_energy" &&
        getParam<Real>("temperature") != defect_temperature)
      paramError("temperature",
                 "must equal the [DefectDynamics] temperature (",
                 defect_temperature,
                 " K) when the solute is coupled to its vacancies, but it is ",
                 getParam<Real>("temperature"),
                 " K.");
    params.set<std::vector<VariableName>>("vacancy_coupled_var") = {
        DefectDynamicsAction::variableName(*vacancy)};
    Real A, C, formation_energy_eV;
    DefectDynamicsAction::resolveCapillaryCoefficients(*vacancy, A, C, formation_energy_eV);
    params.set<Real>("vacancy_formation_energy_eV") = formation_energy_eV;
  }

  _problem->addNodalKernel(kernel_type, name() + "_dyn", params);
}

void
SoluteClusterDynamicsAction::addContinuumPostprocessors()
{
  const std::string var_name = variableName();
  const unsigned int num_cluster_sizes = getParam<unsigned int>("num_cluster_sizes");
  const unsigned int switch_size = getParam<unsigned int>("switch_size");
  const Real initial_width = getParam<Real>("continuum_initial_width");
  const Real width_ratio = getParam<Real>("continuum_width_ratio");

  {
    // The monomer is array component 0 in either representation
    InputParameters params = _factory.getValidParams("ClusterSizeConcentration");
    params.set<std::vector<VariableName>>("clusters") = {var_name};
    params.set<unsigned int>("n_size") = 1;
    _problem->addPostprocessor("ClusterSizeConcentration", name() + "_monomer_concentration", params);
  }

  const Real atomic_volume = getParam<Real>("atomic_volume");
  const Real r1 =
      atomic_volume > 0.0 ? std::cbrt(3.0 * atomic_volume / (4.0 * libMesh::pi)) : 1.0;

  auto addMoment = [&](const std::string & postprocessor_name,
                       const std::string & quantity,
                       unsigned int n_minimum)
  {
    InputParameters params = _factory.getValidParams("ClusterContinuumMoment");
    params.set<std::vector<VariableName>>("clusters") = {var_name};
    params.set<MooseEnum>("quantity") = quantity;
    params.set<unsigned int>("num_cluster_sizes") = num_cluster_sizes;
    params.set<unsigned int>("switch_size") = switch_size;
    params.set<Real>("continuum_initial_width") = initial_width;
    params.set<Real>("continuum_width_ratio") = width_ratio;
    params.set<unsigned int>("n_minimum") = n_minimum;
    if (atomic_volume > 0.0)
      params.set<Real>("atomic_volume") = atomic_volume;
    params.set<Real>("r1") = r1;
    params.set<Real>("radius_exponent") = 1.0 / 3.0;
    _problem->addPostprocessor("ClusterContinuumMoment", postprocessor_name, params);
  };

  addMoment(name() + "_total_density", "density", getParam<unsigned int>("density_n_minimum"));

  const Real resolvable_radius = getParam<Real>("resolvable_radius");
  if (resolvable_radius <= 0.0)
    return;
  if (atomic_volume <= 0.0)
    paramError("resolvable_radius", "requires atomic_volume > 0, to convert the radius to a size.");
  const unsigned int n_minimum =
      static_cast<unsigned int>(std::ceil(std::pow(resolvable_radius / r1, 3.0)));
  addMoment(name() + "_resolvable_density", "density", n_minimum);
  addMoment(name() + "_resolvable_avg_radius", "average_radius", n_minimum);
}

void
SoluteClusterDynamicsAction::actAddPostprocessors()
{
  const std::string var_name = variableName();
  const unsigned int num_cluster_sizes = getParam<unsigned int>("num_cluster_sizes");

  if (usesContinuumTail())
  {
    addContinuumPostprocessors();
    return;
  }
  const bool use_grouping = getParam<bool>("use_grouping");
  const unsigned int explicit_max = getParam<unsigned int>("explicit_max");
  const auto group_spacing = getParam<MooseEnum>("group_spacing");
  const unsigned int group_linear_width = getParam<unsigned int>("group_linear_width");
  const Real group_geometric_ratio = getParam<Real>("group_geometric_ratio");
  const bool enforce_nonnegative = getParam<bool>("enforce_group_nonnegative");
  const Real tolerance_factor = getParam<Real>("group_nonnegative_tolerance_factor");

  {
    InputParameters params = _factory.getValidParams("ClusterSizeConcentrationGrouped");
    params.set<std::vector<VariableName>>("clusters") = {var_name};
    params.set<unsigned int>("n_size") = 1;
    params.set<unsigned int>("num_cluster_sizes") = num_cluster_sizes;
    params.set<bool>("use_grouping") = use_grouping;
    params.set<unsigned int>("explicit_max") = explicit_max;
    params.set<MooseEnum>("group_spacing") = group_spacing;
    params.set<unsigned int>("group_linear_width") = group_linear_width;
    params.set<Real>("group_geometric_ratio") = group_geometric_ratio;
    params.set<bool>("enforce_group_nonnegative") = enforce_nonnegative;
    params.set<Real>("group_nonnegative_tolerance_factor") = tolerance_factor;
    _problem->addPostprocessor(
        "ClusterSizeConcentrationGrouped", name() + "_monomer_concentration", params);
  }
  {
    InputParameters params = _factory.getValidParams("ClusterTotalDensityGrouped");
    params.set<std::vector<VariableName>>("clusters") = {var_name};
    params.set<unsigned int>("n_minimum") = getParam<unsigned int>("density_n_minimum");
    params.set<unsigned int>("num_cluster_sizes") = num_cluster_sizes;
    params.set<bool>("use_grouping") = use_grouping;
    params.set<unsigned int>("explicit_max") = explicit_max;
    params.set<MooseEnum>("group_spacing") = group_spacing;
    params.set<unsigned int>("group_linear_width") = group_linear_width;
    params.set<Real>("group_geometric_ratio") = group_geometric_ratio;
    params.set<bool>("enforce_group_nonnegative") = enforce_nonnegative;
    params.set<Real>("group_nonnegative_tolerance_factor") = tolerance_factor;
    _problem->addPostprocessor("ClusterTotalDensityGrouped", name() + "_total_density", params);
  }

  // Resolvable-cluster density and mean radius, with the same names and meaning as on a
  // [DefectDynamics] species.
  const Real resolvable_radius = getParam<Real>("resolvable_radius");
  if (resolvable_radius <= 0.0)
    return;
  const Real atomic_volume = getParam<Real>("atomic_volume");
  if (atomic_volume <= 0.0)
    paramError("resolvable_radius", "requires atomic_volume > 0, to convert the radius to a size.");

  // Monomer radius of the spherical law the solute kernel uses: r_n = r1 * n^(1/3).
  const Real r1 = std::cbrt(3.0 * atomic_volume / (4.0 * libMesh::pi));
  const unsigned int n_minimum =
      static_cast<unsigned int>(std::ceil(std::pow(resolvable_radius / r1, 3.0)));

  {
    InputParameters params = _factory.getValidParams("ClusterTotalDensityGrouped");
    params.set<std::vector<VariableName>>("clusters") = {var_name};
    params.set<unsigned int>("n_minimum") = n_minimum;
    params.set<Real>("atomic_volume") = atomic_volume;
    params.set<unsigned int>("num_cluster_sizes") = num_cluster_sizes;
    params.set<bool>("use_grouping") = use_grouping;
    params.set<unsigned int>("explicit_max") = explicit_max;
    params.set<MooseEnum>("group_spacing") = group_spacing;
    params.set<unsigned int>("group_linear_width") = group_linear_width;
    params.set<Real>("group_geometric_ratio") = group_geometric_ratio;
    params.set<bool>("enforce_group_nonnegative") = enforce_nonnegative;
    params.set<Real>("group_nonnegative_tolerance_factor") = tolerance_factor;
    _problem->addPostprocessor(
        "ClusterTotalDensityGrouped", name() + "_resolvable_density", params);
  }
  {
    InputParameters params = _factory.getValidParams("ClusterAverageRadiusGrouped");
    params.set<std::vector<VariableName>>("clusters") = {var_name};
    params.set<Real>("r1") = r1;
    params.set<Real>("radius_exponent") = 1.0 / 3.0;
    params.set<unsigned int>("n_minimum") = n_minimum;
    params.set<unsigned int>("num_cluster_sizes") = num_cluster_sizes;
    params.set<bool>("use_grouping") = use_grouping;
    params.set<unsigned int>("explicit_max") = explicit_max;
    params.set<MooseEnum>("group_spacing") = group_spacing;
    params.set<unsigned int>("group_linear_width") = group_linear_width;
    params.set<Real>("group_geometric_ratio") = group_geometric_ratio;
    params.set<bool>("enforce_group_nonnegative") = enforce_nonnegative;
    params.set<Real>("group_nonnegative_tolerance_factor") = tolerance_factor;
    _problem->addPostprocessor(
        "ClusterAverageRadiusGrouped", name() + "_resolvable_avg_radius", params);
  }
}
