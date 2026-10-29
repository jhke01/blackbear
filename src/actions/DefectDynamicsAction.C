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

#include "DefectDynamicsAction.h"
#include "DefectDynamicsGridAction.h"
#include "DefectRecoveryAction.h"
#include "DefectSpeciesAction.h"

#include "FEProblemBase.h"

namespace
{
// Matches the constants hardcoded in GenericVacancyVoidGroupedNodalKernel.C's and
// GenericInterstitialLoopGroupedNodalKernel.C's own anonymous namespaces, so every quantity this
// action derives from them is bit-identical to the value those kernels themselves use, not merely
// numerically close.
constexpr Real defect_dynamics_kB = 1.380649e-23;
constexpr Real defect_dynamics_eV = 1.602176634e-19;

/**
 * Z_loop that makes a loop's size-dependent bias Z(n) = Z + [sqrt(b/(8 pi a0)) Z_loop - Z]
 * n^(-alpha/2) equal to the constant Z at every size: sqrt(b/(8 pi a0)) Z_loop = Z, so Z_loop = Z /
 * sqrt(b/(8 pi a0)). (About Z/0.152 for b = a0/sqrt(3).) The loop's Burgers vector is its own
 * burgers_vector, or a0/sqrt(3) when that is left at zero, as in the kernels.
 */
Real
constantBiasZLoop(Real z, Real burgers_vector_param, Real lattice_parameter)
{
  const Real b =
      burgers_vector_param != 0.0 ? burgers_vector_param : lattice_parameter / std::sqrt(3.0);
  return z / std::sqrt(b / (8.0 * libMesh::pi * lattice_parameter));
}
} // namespace

registerMooseAction("BlackBearApp", DefectDynamicsAction, "add_variable");
registerMooseAction("BlackBearApp", DefectDynamicsAction, "add_function");
registerMooseAction("BlackBearApp", DefectDynamicsAction, "add_ic");
registerMooseAction("BlackBearApp", DefectDynamicsAction, "add_nodal_kernel");
registerMooseAction("BlackBearApp", DefectDynamicsAction, "add_postprocessor");
registerMooseAction("BlackBearApp", DefectDynamicsAction, "add_output");
registerMooseAction("BlackBearApp", DefectDynamicsAction, "add_vector_postprocessor");

InputParameters
DefectDynamicsAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription(
      "A species-organized front end to the grouped defect-dynamics kernels: one or two named "
      "species sub-blocks (e.g. vacancy, interstitial) plus one Grid sub-block build the "
      "grouped array variables, initial conditions, self/opposite-absorption kernels, "
      "recombination kernels, and monomer/resolvable-radius postprocessors, deriving every "
      "coupling quantity the low-level kernels need instead of asking the user to type it.");

  params.addRequiredParam<Real>("temperature", "Temperature T [K], shared by every species.");
  params.addRequiredParam<Real>(
      "atomic_volume", "Atomic volume V_at [m^3], shared by every species.");
  params.addRequiredParam<Real>(
      "lattice_parameter",
      "Lattice parameter a0 [m], shared by every species (also the default source for a loop "
      "species' Burgers vector, b = a0/sqrt(3)).");

  params.addRequiredParam<Real>(
      "dpa_rate", "Displacement-per-atom production rate [dpa/s]. Together with "
      "cascade_efficiency, gives each species' monomer generation rate, "
      "generation = dpa_rate * cascade_efficiency.");
  params.addRequiredParam<Real>(
      "cascade_efficiency",
      "Fraction of dpa_rate that survives cascade recombination to produce free monomers "
      "(a monomer-only cascade spectrum; see [theory/defect_dynamics.md]).");
  params.addParam<Real>(
      "dislocation_density",
      0.0,
      "Dislocation density rho, in m^-2, for the background sink, sink = D*Z*rho, applied to every "
      "species' monomer. Defaults to zero, which turns the constant sink off. With a `[Recovery]` "
      "sub-block the sink instead comes entirely from the time-evolving dislocation density that "
      "block describes (Pokor et al. 2004, Eqs 7-8), and a nonzero dislocation_density is "
      "ignored with a warning.");

  params.addParam<Real>(
      "recombination_radius",
      0.0,
      "Recombination radius r_rec [m] for direct vacancy-interstitial recombination, "
      "recombination_rate = 4*pi*r_rec*(D_i + D_v)/V_at. Only meaningful when both a void and a "
      "loop species are present; if left at zero, no PointDefectRecombinationNodalKernel is "
      "added (opposite-species absorption, independently controlled by opposite_absorption "
      "below, is unaffected by this).");

  params.addParam<bool>(
      "opposite_absorption",
      true,
      "Whether an interstitial monomer absorbed onto a void shrinks it, and a vacancy monomer "
      "absorbed onto a loop shrinks it, when both a void and a loop species are present. "
      "Defaults to on, matching every validated production input; set to false to isolate "
      "recombination (via recombination_radius) from opposite-species absorption, or to run "
      "two species with no coupling between them at all.");

  params.addParam<MooseEnum>(
      "tail",
      MooseEnum("grouping fokker_planck", "grouping"),
      "How the large cluster sizes are stored. 'grouping' (default): bins with a linear "
      "reconstruction, set per species by explicit_max and the group_* parameters. 'fokker_planck': "
      "sizes up to switch_size are individual unknowns and larger sizes are finite-volume cells of "
      "a Fokker-Planck continuum (Pokor et al. 2004, Section 4.1), set by switch_size and the "
      "continuum_* parameters. The continuum needs far fewer unknowns for the same accuracy, but "
      "supports only the capillary binding-energy model, and both species must have the same "
      "maximum cluster size.");
  params.addParam<unsigned int>(
      "switch_size",
      100,
      "Largest cluster size that is an individual unknown when tail = fokker_planck; larger sizes "
      "are continuum cells. Must be at least 4 when a species has cascade seeding, and below every "
      "species' maximum cluster size.");
  params.addParam<Real>("continuum_initial_width",
                        1.0,
                        "Width of the first continuum cell when tail = fokker_planck.");
  params.addParam<Real>("continuum_width_ratio",
                        1.04,
                        "Growth factor of the cell width from one cell to the next when tail = "
                        "fokker_planck.");
  params.addParam<bool>("continuum_limiter",
                        true,
                        "Reconstruct the cell values at a cell face with a minmod-limited slope "
                        "(second order) when tail = fokker_planck. If false the cell averages are "
                        "used (first-order upwind).");

  params.addParam<Real>(
      "timestep_relative_change",
      0.0,
      "If positive, adds a postprocessor '<block name>_timestep_limit' (RelativeChangeTimestepLimit) "
      "that limits the relative change of the resolvable densities of the species with "
      "resolvable_radius to this fraction per step (0.1 is a good start), bounded above by "
      "timestep_cap_fraction * t + timestep_cap_minimum. Name it as the "
      "timestep_limiting_postprocessor of the time stepper. Without it, the resolvable densities "
      "are wrong while they rise by many orders of magnitude, for the grouped and the Fokker-Planck "
      "tail alike. Requires resolvable_radius on at least one species.");
  params.addParam<Real>("timestep_cap_fraction",
                        0.01,
                        "Upper limit of the step in proportion to the time, used with "
                        "timestep_relative_change.");
  params.addParam<Real>("timestep_cap_minimum",
                        0.5,
                        "Upper limit of the step at the start [s], used with "
                        "timestep_relative_change.");
  params.addParam<Real>("timestep_value_floor",
                        1.0,
                        "A resolvable density below this value [1/m^3] is measured against it, so "
                        "that a population that has not formed does not limit the step.");
  params.addParam<std::vector<Real>>(
      "size_distribution_times",
      std::vector<Real>(),
      "Simulation times [s] (NOT dose/dpa) at which to write the full per-size cluster "
      "concentration distribution for every species with output_size_distribution = true. "
      "Shared by every species so they stay directly comparable at the same instants. Required "
      "(and must be non-empty) if any species sets output_size_distribution = true.");

  return params;
}

DefectDynamicsAction::DefectDynamicsAction(const InputParameters & params) : Action(params)
{
}

const DefectDynamicsGridAction *
DefectDynamicsAction::grid() const
{
  // The shared [Grid] sub-block is deprecated: each species now holds its own grouping. When it is
  // present it only supplies a value a species did not set itself.
  const auto grids = _awh.getActions<DefectDynamicsGridAction>();
  if (grids.empty())
    return nullptr;
  if (grids.size() > 1)
    mooseError("[DefectDynamics] found ",
               grids.size(),
               " [Grid] sub-blocks; at most one is "
               "expected.");
  return grids[0];
}

const DefectRecoveryAction *
DefectDynamicsAction::recovery() const
{
  const auto recoveries = _awh.getActions<DefectRecoveryAction>();
  if (recoveries.empty())
    return nullptr;
  if (recoveries.size() > 1)
    mooseError("[DefectDynamics] found ", recoveries.size(), " [Recovery] sub-blocks; at most "
              "one is expected.");
  return recoveries[0];
}

const DefectSpeciesAction *
DefectDynamicsAction::networkLoopSpecies(const DefectRecoveryAction & recovery) const
{
  const DefectSpeciesAction * loops = nullptr;
  if (recovery.isParamValid("interstitial_species"))
  {
    const auto wanted = recovery.getParam<std::string>("interstitial_species");
    for (const auto * sp : species())
      if (sp->name() == wanted)
        loops = sp;
    if (!loops)
      mooseError("[Recovery] interstitial_species '", wanted, "' is not a species sub-block.");
  }
  else
    for (const auto * sp : species())
      if (sp->getParam<MooseEnum>("geometry") == "loop")
      {
        if (loops)
          mooseError("[Recovery] model = cellular_network found two species with geometry = "
                     "loop; name the interstitial one with interstitial_species.");
        loops = sp;
      }
  if (!loops)
    mooseError("[Recovery] model = cellular_network needs a species with geometry = loop, whose "
               "monomer concentration drives the network.");
  if (loops->getParam<MooseEnum>("geometry") != "loop")
    mooseError("[Recovery] interstitial_species must have geometry = loop.");
  return loops;
}

std::vector<const DefectSpeciesAction *>
DefectDynamicsAction::species() const
{
  const auto sp = _awh.getActions<DefectSpeciesAction>();
  if (sp.empty())
    mooseError("[DefectDynamics] requires at least one species sub-block (e.g. [vacancy] or "
              "[interstitial]).");
  if (sp.size() > 2)
    mooseError("[DefectDynamics] found ", sp.size(), " species sub-blocks; at most two "
              "are supported.");
  return sp;
}

bool
DefectDynamicsAction::usesContinuumTail() const
{
  return getParam<MooseEnum>("tail") == "fokker_planck";
}

ClusterContinuumLayout
DefectDynamicsAction::continuumLayout(const DefectSpeciesAction & sp) const
{
  return ClusterContinuumLayout(
      sp.resolvedNumClusterSizes(getParam<Real>("atomic_volume"),
                                 getParam<Real>("lattice_parameter")),
      getParam<unsigned int>("switch_size"),
      getParam<Real>("continuum_initial_width"),
      getParam<Real>("continuum_width_ratio"));
}

void
DefectDynamicsAction::validateTail() const
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

  if (getParam<unsigned int>("switch_size") < 2)
    paramError("switch_size", "must be at least 2.");
  if (getParam<Real>("continuum_initial_width") <= 0.0)
    paramError("continuum_initial_width", "must be positive.");
  if (getParam<Real>("continuum_width_ratio") < 1.0)
    paramError("continuum_width_ratio", "must be at least 1.");

  const Real atomic_volume = getParam<Real>("atomic_volume");
  const Real lattice_parameter = getParam<Real>("lattice_parameter");
  const char * grouping_params[] = {"use_grouping",
                                    "explicit_max",
                                    "group_spacing",
                                    "group_linear_width",
                                    "group_geometric_ratio",
                                    "enforce_group_nonnegative",
                                    "group_nonnegative_tolerance_factor"};
  const auto sp_list = species();
  for (const auto * sp : sp_list)
  {
    if (sp->getParam<MooseEnum>("binding_energy_model") != "capillary")
      mooseError("[DefectDynamics] with tail = fokker_planck supports only the capillary "
                 "binding-energy model, but species '",
                 sp->name(),
                 "' sets binding_energy_model = ",
                 std::string(sp->getParam<MooseEnum>("binding_energy_model")),
                 ".");
    for (const auto * name : grouping_params)
      if (sp->parameters().isParamSetByUser(name))
        mooseError("[DefectDynamics] species '",
                   sp->name(),
                   "' sets ",
                   name,
                   ", which only applies to tail = grouping.");

    const unsigned int n = sp->resolvedNumClusterSizes(atomic_volume, lattice_parameter);
    if (getParam<unsigned int>("switch_size") >= n)
      paramError("switch_size",
                 "must be smaller than the maximum cluster size of every species, but species '",
                 sp->name(),
                 "' has ",
                 n,
                 ".");
    // The resolvable population starts at the size n_min of the smallest resolvable radius. When n_min lies above switch_size it falls
    // inside a continuum cell, which is cut there with the cell average, and the resolvable density is then only as accurate as that cut
    // (a 34 percent error in the void density at switch_size 100 against n_min = 152, 0.3 percent at 160).
    if (sp->getParam<Real>("resolvable_radius") > 0.0)
    {
      Real r1, radius_exponent;
      sp->monomerRadiusLaw(atomic_volume, lattice_parameter, r1, radius_exponent);
      const unsigned int n_min = DefectSpeciesAction::clusterSizeFromRadius(
          sp->getParam<Real>("resolvable_radius"), r1, radius_exponent);
      if (getParam<unsigned int>("switch_size") < n_min)
        paramWarning("switch_size",
                     "is below the size ",
                     n_min,
                     " at which species '",
                     sp->name(),
                     "' becomes resolvable (resolvable_radius), so the cut of its resolvable population falls inside a continuum "
                     "cell and its resolvable density is less accurate. Use switch_size of at least ",
                     n_min,
                     ".");
    }
    const bool cascade = sp->getParam<Real>("cascade_f2") + sp->getParam<Real>("cascade_f3") +
                             sp->getParam<Real>("cascade_f4") >
                         0.0;
    if (cascade && getParam<unsigned int>("switch_size") < 4)
      paramError("switch_size",
                 "must be at least 4 when species '",
                 sp->name(),
                 "' has cascade seeding, so that sizes 2 to 4 are individual unknowns.");
  }
  // The two families see each other's cells, so they need identical layouts
  if (sp_list.size() == 2 &&
      sp_list[0]->resolvedNumClusterSizes(atomic_volume, lattice_parameter) !=
          sp_list[1]->resolvedNumClusterSizes(atomic_volume, lattice_parameter))
    mooseError("[DefectDynamics] with tail = fokker_planck needs the same maximum cluster size for "
               "both species, but '",
               sp_list[0]->name(),
               "' has ",
               sp_list[0]->resolvedNumClusterSizes(atomic_volume, lattice_parameter),
               " and '",
               sp_list[1]->name(),
               "' has ",
               sp_list[1]->resolvedNumClusterSizes(atomic_volume, lattice_parameter),
               ".");
}

Real
DefectDynamicsAction::monomerConcentration(Real Ef_eV) const
{
  const Real temperature = getParam<Real>("temperature");
  return 4.0 * std::exp(-Ef_eV * defect_dynamics_eV / (defect_dynamics_kB * temperature));
}

void
DefectDynamicsAction::capillaryCoefficients(
    Real Ef_eV, Real Eb2_eV, Real exponent, Real & A, Real & C)
{
  A = (Ef_eV - Eb2_eV) / (std::pow(2.0, exponent) - 1.0);
  C = Ef_eV - A;
}

void
DefectDynamicsAction::resolveCapillaryCoefficients(
    const DefectSpeciesAction & sp, Real & A, Real & C, Real & formation_energy_eV)
{
  if (sp.usesDirectCapillaryCoefficients())
  {
    A = sp.getParam<Real>("capillary_coefficient_eV");
    C = sp.getParam<Real>("capillary_constant_eV");
    formation_energy_eV = A + C;
  }
  else
  {
    capillaryCoefficients(
        sp.getParam<Real>("Ef"), sp.getParam<Real>("Eb2"), sp.getParam<Real>("exponent"), A, C);
    formation_energy_eV = sp.getParam<Real>("Ef");
  }
}

std::string
DefectDynamicsAction::variableName(const DefectSpeciesAction & sp)
{
  return sp.name() + "_clusters";
}

std::string
DefectDynamicsAction::functionName(const DefectSpeciesAction & sp)
{
  return sp.name() + "_index_fn";
}

bool
DefectDynamicsAction::anySpeciesWantsSizeDistribution() const
{
  bool any = false;
  for (const auto * sp : species())
    if (sp->getParam<bool>("output_size_distribution"))
      any = true;

  if (any && getParam<std::vector<Real>>("size_distribution_times").empty())
    mooseError("[DefectDynamics] has a species with output_size_distribution = true but "
               "size_distribution_times is empty; set size_distribution_times to the simulation "
               "times [s] at which to write the size distribution.");

  return any;
}

std::string
DefectDynamicsAction::sizeDistributionOutputName() const
{
  return name() + "_size_distribution_csv";
}

void
DefectDynamicsAction::act()
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
  else if (_current_task == "add_output")
    actAddOutputs();
  else if (_current_task == "add_vector_postprocessor")
    actAddVectorPostprocessors();
}

void
DefectDynamicsAction::actAddVariables()
{
  const Real atomic_volume = getParam<Real>("atomic_volume");
  const Real lattice_parameter = getParam<Real>("lattice_parameter");
  const auto grid_action = grid();
  for (const auto * sp : species())
  {
    const unsigned int components =
        usesContinuumTail()
            ? continuumLayout(*sp).componentCount()
            : sp->layout(sp->resolvedNumClusterSizes(atomic_volume, lattice_parameter), grid_action)
                  .componentCount();
    InputParameters params = _factory.getValidParams("ArrayMooseVariable");
    params.set<MooseEnum>("family") = "LAGRANGE";
    params.set<MooseEnum>("order") = "FIRST";
    params.set<unsigned int>("components") = components;
    _problem->addVariable("ArrayMooseVariable", variableName(*sp), params);
  }
}

void
DefectDynamicsAction::actAddFunctions()
{
  for (const auto * sp : species())
  {
    Real A, C, formation_energy_eV;
    resolveCapillaryCoefficients(*sp, A, C, formation_energy_eV);
    const Real ic = monomerConcentration(formation_energy_eV);

    InputParameters params = _factory.getValidParams("PiecewiseConstant");
    params.set<std::vector<Real>>("x") = {0.0, 0.9};
    params.set<std::vector<Real>>("y") = {ic, 0.0};
    params.set<MooseEnum>("direction") = "left";
    _problem->addFunction("PiecewiseConstant", functionName(*sp), params);
  }
}

void
DefectDynamicsAction::actAddICs()
{
  const Real atomic_volume = getParam<Real>("atomic_volume");
  const Real lattice_parameter = getParam<Real>("lattice_parameter");
  const auto grid_action = grid();
  for (const auto * sp : species())
  {
    if (usesContinuumTail())
    {
      // Only the monomer is populated: component 0 reads the function at x = 0 and every other
      // component reads it at x >= 1, where it is zero
      InputParameters params = _factory.getValidParams("ArrayFunctionIndexIC");
      params.set<VariableName>("variable") = variableName(*sp);
      params.set<FunctionName>("function") = functionName(*sp);
      _problem->addInitialCondition("ArrayFunctionIndexIC", sp->name() + "_ic", params);
      continue;
    }
    InputParameters params = _factory.getValidParams("GroupedArrayFunctionIndexIC");
    params.set<VariableName>("variable") = variableName(*sp);
    params.set<FunctionName>("function") = functionName(*sp);
    params.set<unsigned int>("num_cluster_sizes") =
        sp->resolvedNumClusterSizes(atomic_volume, lattice_parameter);
    const auto g = sp->grouping(grid_action);
    params.set<bool>("use_grouping") = g.use_grouping;
    params.set<unsigned int>("explicit_max") = g.explicit_max;
    params.set<MooseEnum>("group_spacing") = g.spacing;
    params.set<unsigned int>("group_linear_width") = g.linear_width;
    params.set<Real>("group_geometric_ratio") = g.geometric_ratio;
    _problem->addInitialCondition("GroupedArrayFunctionIndexIC", sp->name() + "_ic", params);
  }
}

void
DefectDynamicsAction::addContinuumKernel(const DefectSpeciesAction & sp,
                                         const DefectSpeciesAction * opposite,
                                         Real generation,
                                         Real sink,
                                         Real formation_energy_eV,
                                         Real energy_coefficient_eV,
                                         Real energy_constant_eV,
                                         bool cascade)
{
  const Real atomic_volume = getParam<Real>("atomic_volume");
  const Real lattice_parameter = getParam<Real>("lattice_parameter");
  const Real Z = sp.getParam<Real>("Z");
  const bool constant_bias = sp.getParam<MooseEnum>("bias") == "constant";

  InputParameters params = _factory.getValidParams("ContinuumLoopNodalKernel");
  params.set<NonlinearVariableName>("variable") = variableName(sp);

  // With cascade seeding the monomer keeps the rest of the generation rate; the seeded clusters are
  // added by this species' CascadeSpectrumNodalKernel
  Real monomer_generation = generation;
  if (cascade)
    monomer_generation *= 1.0 - sp.getParam<Real>("cascade_f2") - sp.getParam<Real>("cascade_f3") -
                          sp.getParam<Real>("cascade_f4");
  params.set<Real>("generation") = monomer_generation;
  params.set<Real>("sink") = sink;

  params.set<Real>("temperature") = getParam<Real>("temperature");
  params.set<Real>("atomic_volume") = atomic_volume;
  params.set<Real>("lattice_parameter") = lattice_parameter;
  params.set<MooseEnum>("geometry") = sp.getParam<MooseEnum>("geometry");
  params.set<Real>("burgers_vector") = sp.getParam<Real>("burgers_vector");
  params.set<Real>("D0") = sp.getParam<Real>("D0");
  params.set<Real>("Em_eV") = sp.getParam<Real>("Em");
  params.set<Real>("Z") = Z;
  params.set<Real>("Z_loop") =
      constant_bias ? constantBiasZLoop(Z, sp.getParam<Real>("burgers_vector"), lattice_parameter)
                    : sp.getParam<Real>("Z_loop");
  params.set<Real>("alpha") = sp.getParam<Real>("alpha");
  params.set<Real>("formation_energy_eV") = formation_energy_eV;
  params.set<Real>("energy_exponent") = sp.getParam<Real>("exponent");
  params.set<Real>("energy_coefficient_eV") = energy_coefficient_eV;
  params.set<Real>("energy_constant_eV") = energy_constant_eV;

  const auto layout = continuumLayout(sp);
  params.set<unsigned int>("num_cluster_sizes") = layout.numClusterSizes();
  params.set<unsigned int>("switch_size") = layout.switchSize();
  params.set<Real>("continuum_initial_width") = getParam<Real>("continuum_initial_width");
  params.set<Real>("continuum_width_ratio") = getParam<Real>("continuum_width_ratio");
  params.set<bool>("use_limiter") = getParam<bool>("continuum_limiter");

  if (opposite)
  {
    // The same derivation as for the grouped loop kernel: this loop absorbs the opposite monomer
    // with the opposite species' Z as the large-loop limit and this species' own Z_loop_opposite
    // and alpha_opposite; the opposite loops absorb this monomer with this species' Z as the limit
    // and their own Z_loop_opposite and alpha_opposite.
    const Real opposite_Z = opposite->getParam<Real>("Z");
    const Real b = sp.getParam<Real>("burgers_vector");
    const Real opposite_b = opposite->getParam<Real>("burgers_vector");
    const Real default_b = lattice_parameter / std::sqrt(3.0);

    params.set<std::vector<VariableName>>("opposite_coupled_var") = {variableName(*opposite)};
    params.set<MooseEnum>("opposite_geometry") = opposite->getParam<MooseEnum>("geometry");
    params.set<Real>("opposite_D0") = opposite->getParam<Real>("D0");
    params.set<Real>("opposite_Em_eV") = opposite->getParam<Real>("Em");
    params.set<Real>("opposite_Z") = opposite_Z;
    params.set<Real>("opposite_Z_loop") =
        constant_bias ? constantBiasZLoop(opposite_Z, b, lattice_parameter)
                      : sp.getParam<Real>("Z_loop_opposite");
    params.set<Real>("opposite_alpha") = sp.getParam<Real>("alpha_opposite");
    params.set<Real>("opposite_burgers_vector") = opposite_b != 0.0 ? opposite_b : default_b;
    params.set<Real>("partner_Z_loop") =
        opposite->getParam<MooseEnum>("bias") == "constant"
            ? constantBiasZLoop(Z, opposite_b, lattice_parameter)
            : opposite->getParam<Real>("Z_loop_opposite");
    params.set<Real>("partner_alpha") = opposite->getParam<Real>("alpha_opposite");
  }

  _problem->addNodalKernel("ContinuumLoopNodalKernel", sp.name() + "_dyn", params);
}

void
DefectDynamicsAction::actAddNodalKernels()
{
  const auto grid_action = grid();
  const auto recovery_action = recovery();
  const auto sp_list = species();
  const Real temperature = getParam<Real>("temperature");
  const Real atomic_volume = getParam<Real>("atomic_volume");
  const Real lattice_parameter = getParam<Real>("lattice_parameter");
  const Real dislocation_density = getParam<Real>("dislocation_density");
  const Real generation = getParam<Real>("dpa_rate") * getParam<Real>("cascade_efficiency");
  if (getParam<Real>("temperature") <= 0.0)
    paramError("temperature", "must be positive.");
  if (getParam<Real>("dpa_rate") < 0.0)
    paramError("dpa_rate", "must not be negative.");
  if (getParam<Real>("cascade_efficiency") < 0.0 || getParam<Real>("cascade_efficiency") > 1.0)
    paramError("cascade_efficiency", "must be between 0 and 1.");
  if (dislocation_density < 0.0)
    paramError("dislocation_density", "must not be negative; use 0 to turn the sink off.");
  if (getParam<Real>("recombination_radius") < 0.0)
    paramError("recombination_radius", "must not be negative; use 0 to turn recombination off.");

  // A species activates cascade seeding (Pokor et al. 2004, Eqs 1a-1e) simply by setting a
  // nonzero cascade_f2/cascade_f3/cascade_f4 -- the same "nonzero means on" idiom already used by
  // recombination_radius and dislocation_density, independent per species.
  auto cascadeActive = [](const DefectSpeciesAction & sp)
  {
    return sp.getParam<Real>("cascade_f2") + sp.getParam<Real>("cascade_f3") +
               sp.getParam<Real>("cascade_f4") >
           0.0;
  };

  // Cascade seeding (sizes 2-4) needs explicit (ungrouped) components for those sizes to exist
  // at all; checked per species (its own grouping and resolved cap) for whichever species has
  // cascade seeding active.
  for (const auto * sp : sp_list)
  {
    if (!cascadeActive(*sp))
      continue;
    if (usesContinuumTail())
      continue; // switch_size >= 4 was checked in validateTail()
    const auto g = sp->grouping(grid_action);
    if (g.use_grouping && g.explicit_max < 4)
      mooseError("[DefectDynamics] requires explicit_max >= 4 when species '",
                 sp->name(),
                 "' has cascade seeding active, but explicit_max = ",
                 g.explicit_max,
                 ".");
    if (sp->resolvedNumClusterSizes(atomic_volume, lattice_parameter) < 4)
      mooseError("[DefectDynamics] requires species '", sp->name(), "' to represent at least "
                "cluster size 4 when its own cascade_f2/cascade_f3/cascade_f4 are active.");
  }

  // [Recovery], when present, replaces the constant dislocation sink with a time-evolving one
  // (Pokor et al. 2004, Eqs 7-8 use one un-subscripted dislocation density for both species'
  // sink terms) -- it does not add to dislocation_density, since Pokor's model has exactly one
  // sink term per species, not two. A nonzero dislocation_density alongside [Recovery] would
  // otherwise be silently unused, so warn explicitly instead. (dislocation_density = 0, with or
  // without [Recovery], remains the existing, documented way to turn the constant sink off.)
  if (recovery_action && dislocation_density > 0.0)
    paramWarning("dislocation_density",
                "dislocation_density is unused when a [Recovery] sub-block is given: the "
                "entire dislocation sink comes from [Recovery]'s own rho0 instead.");

  // Exactly one void species is still required when two species are present: Pokor's own
  // dislocation-bias absorption model (Eqs 4a-4f) treats BOTH clusters identically (same loop
  // form), so two loop species already has a well-defined, validated meaning (vacancy modeled as
  // a loop, matching Pokor et al. 2004's own equations -- see
  // sa304_osiris_vacancy_as_loop.i and defect_dynamics_generalization_plan.md Section 9.9/4a).
  // Two VOID species has no such precedent or physical motivation under this codebase's models
  // (finite-sink-strength absorption, Eqs 2a-2f, is specifically a 3-D vacancy-void law) and
  // remains rejected.
  unsigned int void_count = 0;
  for (const auto * sp : sp_list)
    if (sp->getParam<MooseEnum>("geometry") == "void")
      ++void_count;
  if (sp_list.size() == 2 && void_count == 2)
    mooseError("[DefectDynamics] supports at most one void species; found two void species. Two "
              "loop species (e.g. vacancy modeled as a loop, following Pokor et al. 2004's own "
              "equations) are supported; two void species are not.");

  // "The other species" -- used throughout below for opposite-species coupling and
  // recombination. Only ever called when exactly two species are present.
  auto other = [&](const DefectSpeciesAction & sp) -> const DefectSpeciesAction *
  {
    for (const auto * candidate : sp_list)
      if (candidate != &sp)
        return candidate;
    return nullptr;
  };

  // Both species present is necessary but not sufficient for either coupling term: opposite-
  // species absorption is gated on opposite_absorption (on by default, matching every validated
  // production input), and point-defect recombination is a SEPARATE, independently opt-in
  // physical process gated on recombination_radius being given -- a user may legitimately want
  // one without the other (e.g. to isolate one physical effect at a time, as
  // defect_dynamics_opposite_absorption_small.i does for recombination, or to turn off
  // absorption while keeping recombination).
  const bool coupled = sp_list.size() == 2;
  const bool opposite_absorption_active = coupled && getParam<bool>("opposite_absorption");
  const bool recombination_active = coupled && getParam<Real>("recombination_radius") > 0.0;

  // D*exp(-Em/kT) for each species present, needed both for this species' own sink and (when
  // recombination is active) for the recombination rate.
  auto diffusivity = [&](const DefectSpeciesAction & sp)
  {
    return sp.getParam<Real>("D0") *
           std::exp(-sp.getParam<Real>("Em") * defect_dynamics_eV / (defect_dynamics_kB * temperature));
  };

  Real recombination_rate = 0.0;
  if (recombination_active)
    recombination_rate = 4.0 * libMesh::pi * getParam<Real>("recombination_radius") *
                        (diffusivity(*sp_list[0]) + diffusivity(*sp_list[1])) / atomic_volume;

  for (const auto * sp : sp_list)
  {
    const std::string var_name = variableName(*sp);

    // Time derivative, shared shape for every species.
    {
      InputParameters params = _factory.getValidParams("ArrayTimeDerivativeNodalKernel");
      params.set<NonlinearVariableName>("variable") = var_name;
      _problem->addNodalKernel("ArrayTimeDerivativeNodalKernel", sp->name() + "_dot", params);
    }

    const bool is_void = sp->getParam<MooseEnum>("geometry") == "void";
    const DefectSpeciesAction * opposite =
        opposite_absorption_active ? other(*sp) : nullptr;

    // Two mutually exclusive ways to get the capillary law's A/C: derived from a monomer
    // formation energy and a dimer binding energy (the default, physically-intuitive path, used
    // by every validated production input), or given directly (for a user who already has a
    // fitted/literature power law, e.g. Soneda's A*n^(2/3) with C = 0). Either way, both paths
    // converge on the same three kernel inputs: A, C, and formation_energy_eV = A + C, which is
    // E_f(n) evaluated at n = 1 -- exactly what Ef names in the other path.
    Real A, C, formation_energy_eV;
    resolveCapillaryCoefficients(*sp, A, C, formation_energy_eV);

    if (usesContinuumTail())
    {
      addContinuumKernel(*sp,
                         opposite,
                         generation,
                         recovery_action ? 0.0 : diffusivity(*sp) * sp->getParam<Real>("Z") * dislocation_density,
                         formation_energy_eV,
                         A,
                         C,
                         cascadeActive(*sp));
      continue;
    }

    const std::string kernel_type =
        is_void ? "VacancyVoidGroupedNodalKernel" : "InterstitialLoopGroupedNodalKernel";
    InputParameters params = _factory.getValidParams(kernel_type);
    params.set<NonlinearVariableName>("variable") = var_name;

    params.set<unsigned int>("num_cluster_sizes") =
        sp->resolvedNumClusterSizes(atomic_volume, lattice_parameter);
    // This species' OWN array grouping (its variable, its kernel and its IC all use it).
    const auto own_grouping = sp->grouping(grid_action);
    params.set<bool>("use_grouping") = own_grouping.use_grouping;
    params.set<unsigned int>("explicit_max") = own_grouping.explicit_max;
    params.set<MooseEnum>("group_spacing") = own_grouping.spacing;
    params.set<unsigned int>("group_linear_width") = own_grouping.linear_width;
    params.set<Real>("group_geometric_ratio") = own_grouping.geometric_ratio;
    params.set<bool>("enforce_group_nonnegative") = own_grouping.enforce_nonnegative;
    params.set<Real>("group_nonnegative_tolerance_factor") = own_grouping.tolerance_factor;

    // When this species has cascade seeding active, its own CascadeSpectrumNodalKernel instance
    // (added below) seeds cluster sizes 2-4 directly; the monomer's own share of the generation
    // rate shrinks by the same fractions so total point-defect production stays conserved at
    // exactly "generation" (Pokor et al. 2004, Eqs 1a-1e).
    Real monomer_generation = generation;
    if (cascadeActive(*sp))
      monomer_generation *= 1.0 - sp->getParam<Real>("cascade_f2") -
                            sp->getParam<Real>("cascade_f3") - sp->getParam<Real>("cascade_f4");
    params.set<Real>("generation") = monomer_generation;
    // When [Recovery] is present, this species' DislocationRecoveryNodalKernel instance (added
    // below) supplies the entire time-evolving dislocation sink instead; the kernel's own
    // constant "sink" param is left at its default of zero rather than also being set here.
    if (!recovery_action)
      params.set<Real>("sink") = diffusivity(*sp) * sp->getParam<Real>("Z") * dislocation_density;

    params.set<Real>("temperature") = temperature;
    params.set<Real>("atomic_volume") = atomic_volume;
    params.set<Real>("lattice_parameter") = lattice_parameter;
    params.set<Real>("formation_energy_eV") = formation_energy_eV;

    params.set<Real>(is_void ? "D0_v" : "D0_i") = sp->getParam<Real>("D0");
    params.set<Real>(is_void ? "Em_v_eV" : "Em_i_eV") = sp->getParam<Real>("Em");
    params.set<Real>(is_void ? "void_energy_exponent" : "loop_energy_exponent") =
        sp->getParam<Real>("exponent");
    params.set<Real>(is_void ? "void_energy_coefficient_eV" : "loop_energy_coefficient_eV") = A;
    params.set<Real>(is_void ? "void_energy_constant_eV" : "loop_energy_constant_eV") = C;
    params.set<MooseEnum>("binding_energy_model") = sp->getParam<MooseEnum>("binding_energy_model");
    params.set<std::vector<Real>>("binding_energy_table_eV") =
        sp->getParam<std::vector<Real>>("binding_energy_table_eV");
    params.set<std::vector<Real>>("formation_energy_table_eV") =
        sp->getParam<std::vector<Real>>("formation_energy_table_eV");

    if (!is_void)
    {
      const MooseEnum bias = sp->getParam<MooseEnum>("bias");
      const Real Z = sp->getParam<Real>("Z");
      params.set<Real>("z_i") = Z;
      params.set<Real>("z_loop_i") =
          bias == "constant"
              ? constantBiasZLoop(Z, sp->getParam<Real>("burgers_vector"), lattice_parameter)
              : sp->getParam<Real>("Z_loop");
      params.set<Real>("alpha_loop_i") = sp->getParam<Real>("alpha");
      params.set<Real>("burgers_vector") = sp->getParam<Real>("burgers_vector");
    }

    if (opposite)
    {
      params.set<Real>(is_void ? "D0_i" : "D0_v") = opposite->getParam<Real>("D0");
      params.set<Real>(is_void ? "Em_i_eV" : "Em_v_eV") = opposite->getParam<Real>("Em");

      params.set<std::vector<VariableName>>("opposite_coupled_var") = {variableName(*opposite)};
      params.set<unsigned int>("opposite_num_cluster_sizes") =
          opposite->resolvedNumClusterSizes(atomic_volume, lattice_parameter);
      // The kernel rebuilds the PARTNER's array to sum its absorption (Piece B), so it must be
      // given the partner's own grouping, not this species'.
      const auto partner_grouping = opposite->grouping(grid_action);
      params.set<bool>("opposite_use_grouping") = partner_grouping.use_grouping;
      params.set<unsigned int>("opposite_explicit_max") = partner_grouping.explicit_max;
      params.set<MooseEnum>("opposite_group_spacing") = partner_grouping.spacing;
      params.set<unsigned int>("opposite_group_linear_width") = partner_grouping.linear_width;
      params.set<Real>("opposite_group_geometric_ratio") = partner_grouping.geometric_ratio;
      params.set<bool>("opposite_enforce_group_nonnegative") = partner_grouping.enforce_nonnegative;
      params.set<Real>("opposite_group_nonnegative_tolerance_factor") =
          partner_grouping.tolerance_factor;

      if (is_void)
      {
        // The void's own absorption of the opposite (interstitial) species uses the LOOP's
        // vacancy-absorption bias law (oppositeBias(), the same functional form as the loop's own
        // loopBias()): the n -> infinity limit is the VACANCY species' own Z (this kernel's "sp"),
        // and the loop species' Z_loop_opposite/alpha_opposite are the coefficient/exponent that
        // set how it approaches that limit for a loop absorbing a vacancy monomer specifically
        // (see DefectSpeciesAction).
        params.set<Real>("opposite_z_v") = sp->getParam<Real>("Z");
        params.set<Real>("opposite_z_loop_v") =
            opposite->getParam<MooseEnum>("bias") == "constant"
                ? constantBiasZLoop(sp->getParam<Real>("Z"),
                                    opposite->getParam<Real>("burgers_vector"),
                                    lattice_parameter)
                : opposite->getParam<Real>("Z_loop_opposite");
        params.set<Real>("opposite_alpha_loop_v") = opposite->getParam<Real>("alpha_opposite");
        params.set<Real>("opposite_burgers_vector") = opposite->getParam<Real>("burgers_vector");
      }
      else
      {
        // The loop's own absorption of the opposite (vacancy) species uses the loop's OWN
        // vacancy-absorption bias (z_v/z_loop_v/alpha_loop_v, the same oppositeBias() form): the
        // n -> infinity limit is the VACANCY species' own Z (this kernel's "opposite"), and this
        // loop's own Z_loop_opposite/alpha_opposite are the coefficient/exponent.
        params.set<Real>("z_v") = opposite->getParam<Real>("Z");
        params.set<Real>("z_loop_v") = sp->getParam<MooseEnum>("bias") == "constant"
                                           ? constantBiasZLoop(opposite->getParam<Real>("Z"),
                                                               sp->getParam<Real>("burgers_vector"),
                                                               lattice_parameter)
                                           : sp->getParam<Real>("Z_loop_opposite");
        params.set<Real>("alpha_loop_v") = sp->getParam<Real>("alpha_opposite");

        // When the opposite species is ALSO a loop (e.g. vacancy modeled as a loop, following
        // Pokor et al. 2004 Eqs 4a-4f), the monomer-pool loss (Piece B) must sum the opposite
        // loops' absorption of THIS species' monomer with the loop form, not the void form: it is
        // the very coefficient the opposite kernel evaluates as its own Piece A (its z_v,
        // z_loop_v, alpha_loop_v, set by the branch above on that species' pass), so both sides of
        // one absorption event use one rate. A void opposite keeps the kernel default (void form).
        if (opposite->getParam<MooseEnum>("geometry") == "loop")
        {
          params.set<MooseEnum>("opposite_geometry") = "loop";
          const Real opposite_burgers = opposite->getParam<Real>("burgers_vector");
          params.set<Real>("opposite_burgers_vector") =
              opposite_burgers != 0.0 ? opposite_burgers : lattice_parameter / std::sqrt(3.0);
          params.set<Real>("opposite_z_i") = sp->getParam<Real>("Z");
          params.set<Real>("opposite_z_loop_i") =
              opposite->getParam<MooseEnum>("bias") == "constant"
                  ? constantBiasZLoop(sp->getParam<Real>("Z"), opposite_burgers, lattice_parameter)
                  : opposite->getParam<Real>("Z_loop_opposite");
          params.set<Real>("opposite_alpha_loop_i") = opposite->getParam<Real>("alpha_opposite");
        }
      }
    }

    _problem->addNodalKernel(kernel_type, sp->name() + "_dyn", params);
  }

  if (recovery_action)
  {
    for (const auto * sp : sp_list)
    {
      InputParameters params = _factory.getValidParams("DislocationRecoveryNodalKernel");
      params.set<NonlinearVariableName>("variable") = variableName(*sp);
      params.set<Real>("D0") = sp->getParam<Real>("D0");
      params.set<Real>("Em_eV") = sp->getParam<Real>("Em");
      params.set<Real>("temperature") = temperature;
      params.set<Real>("Z") = sp->getParam<Real>("Z");
      params.set<Real>("burgers_vector") = recovery_action->getParam<Real>("burgers_vector");
      if (recovery_action->cellularNetwork())
        params.set<PostprocessorName>("dislocation_density_postprocessor") =
            name() + "_network_density";
      else
      {
        params.set<Real>("rho0") = recovery_action->getParam<Real>("rho0");
        params.set<Real>("K") = recovery_action->getParam<Real>("K");
      }
      _problem->addNodalKernel(
          "DislocationRecoveryNodalKernel", sp->name() + "_recovery", params);
    }

    if (recovery_action->cellularNetwork() && recovery_action->getParam<bool>("loop_transfer"))
    {
      if (!usesContinuumTail())
        mooseError("[Recovery] loop_transfer = true needs tail = fokker_planck: the line length of "
                   "the loops that join the network is summed over the continuum layout.");
      const auto * loops = networkLoopSpecies(*recovery_action);
      InputParameters params = _factory.getValidParams("LoopNetworkTransferNodalKernel");
      params.set<NonlinearVariableName>("variable") = variableName(*loops);
      params.set<PostprocessorName>("dislocation_density_postprocessor") =
          name() + "_network_density";
      params.set<Real>("D0") = loops->getParam<Real>("D0");
      params.set<Real>("Em_eV") = loops->getParam<Real>("Em");
      params.set<Real>("temperature") = temperature;
      params.set<Real>("burgers_vector") = recovery_action->getParam<Real>("burgers_vector");
      params.set<Real>("bias") = recovery_action->getParam<Real>("bias");
      params.set<Real>("ca_crho") = recovery_action->getParam<Real>("ca_crho");
      params.set<Real>("poisson_ratio") = recovery_action->getParam<Real>("poisson_ratio");
      params.set<unsigned int>("minimum_size") =
          recovery_action->getParam<unsigned int>("transfer_minimum_size");
      _problem->addNodalKernel("LoopNetworkTransferNodalKernel", loops->name() + "_network_transfer", params);
    }
  }

  for (const auto * sp : sp_list)
  {
    if (!cascadeActive(*sp))
      continue;
    InputParameters params = _factory.getValidParams("CascadeSpectrumNodalKernel");
    params.set<NonlinearVariableName>("variable") = variableName(*sp);
    params.set<Real>("generation") = generation;
    params.set<Real>("f2") = sp->getParam<Real>("cascade_f2");
    params.set<Real>("f3") = sp->getParam<Real>("cascade_f3");
    params.set<Real>("f4") = sp->getParam<Real>("cascade_f4");
    _problem->addNodalKernel("CascadeSpectrumNodalKernel", sp->name() + "_cascade", params);
  }

  if (recombination_active)
  {
    for (const auto * sp : sp_list)
    {
      const auto * opp = other(*sp);
      InputParameters params = _factory.getValidParams("PointDefectRecombinationNodalKernel");
      params.set<NonlinearVariableName>("variable") = variableName(*sp);
      params.set<std::vector<VariableName>>("coupled_var") = {variableName(*opp)};
      params.set<Real>("recombination_rate") = recombination_rate;
      _problem->addNodalKernel(
          "PointDefectRecombinationNodalKernel", sp->name() + "_recombination", params);
    }
  }
}

namespace
{
/// Copies the grid shape and nonnegative-reconstruction parameters shared by every grouped
/// postprocessor/vector postprocessor onto params, so each caller only needs to set the
/// object-specific parameters (clusters, num_cluster_sizes, n_minimum/n_size, r1, etc.).
void
setGroupingLayoutParams(InputParameters & params, const DefectSpeciesAction::Grouping & g)
{
  params.set<bool>("use_grouping") = g.use_grouping;
  params.set<unsigned int>("explicit_max") = g.explicit_max;
  params.set<MooseEnum>("group_spacing") = g.spacing;
  params.set<unsigned int>("group_linear_width") = g.linear_width;
  params.set<Real>("group_geometric_ratio") = g.geometric_ratio;
  params.set<bool>("enforce_group_nonnegative") = g.enforce_nonnegative;
  params.set<Real>("group_nonnegative_tolerance_factor") = g.tolerance_factor;
}
} // namespace

void
DefectDynamicsAction::actAddPostprocessors()
{
  const auto grid_action = grid();
  const Real atomic_volume = getParam<Real>("atomic_volume");
  const Real lattice_parameter = getParam<Real>("lattice_parameter");

  const auto recovery_action = recovery();
  const bool cellular = recovery_action && recovery_action->cellularNetwork();
  if (cellular)
  {
    // Rest (2006): the network density is driven by the interstitial loop species' monomers
    const DefectSpeciesAction * loops = networkLoopSpecies(*recovery_action);
    InputParameters params = _factory.getValidParams("CellularNetworkDislocationDensity");
    params.set<PostprocessorName>("interstitial_concentration") =
        loops->name() + "_monomer_concentration";
    params.set<Real>("D0") = loops->getParam<Real>("D0");
    params.set<Real>("Em_eV") = loops->getParam<Real>("Em");
    params.set<Real>("temperature") = getParam<Real>("temperature");
    params.set<Real>("burgers_vector") = recovery_action->getParam<Real>("burgers_vector");
    params.set<Real>("atomic_volume") = atomic_volume;
    params.set<Real>("rho0") = recovery_action->getParam<Real>("rho0");
    params.set<Real>("bias") = recovery_action->getParam<Real>("bias");
    params.set<Real>("ca_crho") = recovery_action->getParam<Real>("ca_crho");
    params.set<Real>("poisson_ratio") = recovery_action->getParam<Real>("poisson_ratio");
    params.set<Real>("grain_size") = recovery_action->getParam<Real>("grain_size");
    if (recovery_action->getParam<bool>("loop_transfer"))
    {
      if (!usesContinuumTail())
        mooseError("[Recovery] loop_transfer = true needs tail = fokker_planck: the line length of "
                   "the loops that join the network is summed over the continuum layout.");
      // The loops that join the network: all loops from transfer_minimum_size up, counted with the
      // same layout as the cluster variable, so that the line length removed by
      // LoopNetworkTransferNodalKernel is exactly the line length added here
      Real r1, radius_exponent;
      loops->monomerRadiusLaw(atomic_volume, getParam<Real>("lattice_parameter"), r1, radius_exponent);
      InputParameters line = _factory.getValidParams("ClusterContinuumMoment");
      line.set<std::vector<VariableName>>("clusters") = {variableName(*loops)};
      line.set<MooseEnum>("quantity") = "line_length";
      line.set<unsigned int>("num_cluster_sizes") =
          loops->resolvedNumClusterSizes(atomic_volume, getParam<Real>("lattice_parameter"));
      line.set<unsigned int>("switch_size") = getParam<unsigned int>("switch_size");
      line.set<Real>("continuum_initial_width") = getParam<Real>("continuum_initial_width");
      line.set<Real>("continuum_width_ratio") = getParam<Real>("continuum_width_ratio");
      line.set<unsigned int>("n_minimum") = recovery_action->getParam<unsigned int>("transfer_minimum_size");
      line.set<Real>("atomic_volume") = atomic_volume;
      line.set<Real>("r1") = r1;
      line.set<Real>("radius_exponent") = radius_exponent;
      line.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
      _problem->addPostprocessor("ClusterContinuumMoment", name() + "_network_line_length", line);
      params.set<PostprocessorName>("line_length_postprocessor") = name() + "_network_line_length";
    }
    _problem->addPostprocessor("CellularNetworkDislocationDensity", name() + "_network_density", params);
  }

  if (getParam<Real>("timestep_relative_change") > 0.0)
  {
    std::vector<PostprocessorName> monitored;
    if (cellular)
      monitored.push_back(name() + "_network_density");
    for (const auto * sp : species())
      if (sp->getParam<Real>("resolvable_radius") > 0.0)
        monitored.push_back(sp->name() + "_resolvable_density");
    if (monitored.empty())
      paramError("timestep_relative_change",
                 "needs resolvable_radius on at least one species, whose resolvable density it "
                 "limits.");
    InputParameters params = _factory.getValidParams("RelativeChangeTimestepLimit");
    params.set<std::vector<PostprocessorName>>("postprocessors") = monitored;
    params.set<Real>("relative_change_target") = getParam<Real>("timestep_relative_change");
    params.set<Real>("value_floor") = getParam<Real>("timestep_value_floor");
    params.set<Real>("time_fraction") = getParam<Real>("timestep_cap_fraction");
    params.set<Real>("minimum_cap") = getParam<Real>("timestep_cap_minimum");
    _problem->addPostprocessor("RelativeChangeTimestepLimit", name() + "_timestep_limit", params);
  }

  for (const auto * sp : species())
  {
    const std::string var_name = variableName(*sp);
    const unsigned int num_cluster_sizes =
        sp->resolvedNumClusterSizes(atomic_volume, lattice_parameter);

    if (usesContinuumTail())
    {
      // The monomer is array component 0 in either representation
      InputParameters params = _factory.getValidParams("ClusterSizeConcentration");
      params.set<std::vector<VariableName>>("clusters") = {var_name};
      params.set<unsigned int>("n_size") = 1;
      _problem->addPostprocessor(
          "ClusterSizeConcentration", sp->name() + "_monomer_concentration", params);
    }
    else
    {
      InputParameters params = _factory.getValidParams("ClusterSizeConcentrationGrouped");
      params.set<std::vector<VariableName>>("clusters") = {var_name};
      params.set<unsigned int>("n_size") = 1;
      params.set<unsigned int>("num_cluster_sizes") = num_cluster_sizes;
      setGroupingLayoutParams(params, sp->grouping(grid_action));
      _problem->addPostprocessor(
          "ClusterSizeConcentrationGrouped", sp->name() + "_monomer_concentration", params);
    }

    const Real resolvable_radius = sp->getParam<Real>("resolvable_radius");
    if (resolvable_radius <= 0.0)
      continue;

    Real r1, radius_exponent;
    sp->monomerRadiusLaw(atomic_volume, lattice_parameter, r1, radius_exponent);
    const unsigned int n_minimum =
        DefectSpeciesAction::clusterSizeFromRadius(resolvable_radius, r1, radius_exponent);

    if (usesContinuumTail())
    {
      const std::pair<const char *, const char *> quantities[] = {
          {"density", "_resolvable_density"}, {"average_radius", "_resolvable_avg_radius"}};
      for (const auto & [quantity, suffix] : quantities)
      {
        InputParameters params = _factory.getValidParams("ClusterContinuumMoment");
        params.set<std::vector<VariableName>>("clusters") = {var_name};
        params.set<MooseEnum>("quantity") = quantity;
        params.set<unsigned int>("num_cluster_sizes") = num_cluster_sizes;
        params.set<unsigned int>("switch_size") = getParam<unsigned int>("switch_size");
        params.set<Real>("continuum_initial_width") = getParam<Real>("continuum_initial_width");
        params.set<Real>("continuum_width_ratio") = getParam<Real>("continuum_width_ratio");
        params.set<unsigned int>("n_minimum") = n_minimum;
        params.set<Real>("atomic_volume") = atomic_volume;
        params.set<Real>("r1") = r1;
        params.set<Real>("radius_exponent") = radius_exponent;
        _problem->addPostprocessor("ClusterContinuumMoment", sp->name() + suffix, params);
      }
      continue;
    }

    {
      InputParameters params = _factory.getValidParams("ClusterTotalDensityGrouped");
      params.set<std::vector<VariableName>>("clusters") = {var_name};
      params.set<unsigned int>("n_minimum") = n_minimum;
      params.set<Real>("atomic_volume") = atomic_volume;
      params.set<unsigned int>("num_cluster_sizes") = num_cluster_sizes;
      setGroupingLayoutParams(params, sp->grouping(grid_action));
      _problem->addPostprocessor(
          "ClusterTotalDensityGrouped", sp->name() + "_resolvable_density", params);
    }
    {
      InputParameters params = _factory.getValidParams("ClusterAverageRadiusGrouped");
      params.set<std::vector<VariableName>>("clusters") = {var_name};
      params.set<Real>("r1") = r1;
      params.set<Real>("radius_exponent") = radius_exponent;
      params.set<unsigned int>("n_minimum") = n_minimum;
      params.set<unsigned int>("num_cluster_sizes") = num_cluster_sizes;
      setGroupingLayoutParams(params, sp->grouping(grid_action));
      _problem->addPostprocessor(
          "ClusterAverageRadiusGrouped", sp->name() + "_resolvable_avg_radius", params);
    }
  }
}

void
DefectDynamicsAction::actAddOutputs()
{
  if (!anySpeciesWantsSizeDistribution())
    return;

  InputParameters params = _factory.getValidParams("CSV");
  params.set<std::vector<Real>>("sync_times") = getParam<std::vector<Real>>("size_distribution_times");
  params.set<bool>("sync_only") = true;
  // Without this, this output falls back to the common [Outputs] block's file_base and collides
  // with the main CSV output (both would try to write <file_base>.csv).
  params.set<bool>("append_object_name") = true;
  // Restrict this dedicated output to vector postprocessors only, so it doesn't ALSO duplicate
  // the regular scalar postprocessor time history the main [Outputs] csv = true block already
  // writes every timestep.
  params.set<ExecFlagEnum>("execute_postprocessors_on") = EXEC_NONE;
  params.set<ExecFlagEnum>("execute_scalars_on") = EXEC_NONE;
  _problem->addOutput("CSV", sizeDistributionOutputName(), params);
}

void
DefectDynamicsAction::actAddVectorPostprocessors()
{
  if (!anySpeciesWantsSizeDistribution())
    return;

  const auto grid_action = grid();
  const Real atomic_volume = getParam<Real>("atomic_volume");
  const Real lattice_parameter = getParam<Real>("lattice_parameter");
  const std::string output_name = sizeDistributionOutputName();

  for (const auto * sp : species())
  {
    if (!sp->getParam<bool>("output_size_distribution"))
      continue;

    const std::string var_name = variableName(*sp);
    const unsigned int num_cluster_sizes =
        sp->resolvedNumClusterSizes(atomic_volume, lattice_parameter);

    Real r1, radius_exponent;
    sp->monomerRadiusLaw(atomic_volume, lattice_parameter, r1, radius_exponent);

    if (usesContinuumTail())
    {
      InputParameters params = _factory.getValidParams("ClusterContinuumDistribution");
      params.set<std::vector<VariableName>>("clusters") = {var_name};
      params.set<unsigned int>("num_cluster_sizes") = num_cluster_sizes;
      params.set<unsigned int>("switch_size") = getParam<unsigned int>("switch_size");
      params.set<Real>("continuum_initial_width") = getParam<Real>("continuum_initial_width");
      params.set<Real>("continuum_width_ratio") = getParam<Real>("continuum_width_ratio");
      params.set<Real>("r1") = r1;
      params.set<Real>("radius_exponent") = radius_exponent;
      params.set<std::vector<OutputName>>("outputs") = {output_name};
      _problem->addVectorPostprocessor(
          "ClusterContinuumDistribution", sp->name() + "_size_distribution", params);
      continue;
    }

    InputParameters params = _factory.getValidParams("ClusterSizeDistributionGrouped");
    params.set<std::vector<VariableName>>("clusters") = {var_name};
    params.set<unsigned int>("num_cluster_sizes") = num_cluster_sizes;
    params.set<Real>("r1") = r1;
    params.set<Real>("radius_exponent") = radius_exponent;
    setGroupingLayoutParams(params, sp->grouping(grid_action));
    params.set<std::vector<OutputName>>("outputs") = {output_name};
    _problem->addVectorPostprocessor(
        "ClusterSizeDistributionGrouped", sp->name() + "_size_distribution", params);
  }
}
