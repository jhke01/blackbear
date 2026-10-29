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

#include "DefectSpeciesAction.h"
#include "DefectDynamicsGridAction.h"

registerMooseAction("BlackBearApp", DefectSpeciesAction, "meta_action");

InputParameters
DefectSpeciesAction::validParams()
{
  InputParameters params = Action::validParams();
  MooseEnum geometry("void loop");
  MooseEnum bias("constant size_dependent", "size_dependent");
  params.addClassDescription(
      "One defect species (e.g. vacancy or interstitial) under the enclosing DefectDynamics "
      "block: its own diffusivity, capillary formation-energy law, and (for a loop) dislocation "
      "bias. DefectDynamicsAction reads this block, together with the other species present, to "
      "build the grouped array variable, initial condition, kernels, and postprocessors for this "
      "species -- every parameter the underlying kernels need about the OTHER species "
      "('opposite_*') is derived from that other species' own block, never typed here.");

  params.addRequiredParam<MooseEnum>(
      "geometry",
      geometry,
      "Cluster geometry: 'void' for a three-dimensional vacancy cluster (radius ~ n^(1/3), "
      "finite-sink-strength absorption correction) or 'loop' for a two-dimensional interstitial "
      "cluster (radius ~ sqrt(n), dislocation-bias absorption).");

  params.addRequiredParam<Real>("D0", "Diffusion prefactor D0 [m^2/s] for this species' monomer.");
  params.addRequiredParam<Real>("Em", "Migration energy Em [eV] for this species' monomer.");

  params.addParam<Real>(
      "Ef",
      0.0,
      "Monomer formation energy E_f [eV]: used (together with Eb2) to derive the capillary "
      "formation-energy law's coefficient and constant, A = (Ef - Eb2)/(2^exponent - 1), "
      "C = Ef - A, and to derive this species' initial monomer concentration, 4*exp(-Ef/kT). "
      "The default, physically-intuitive path: set this and Eb2, and leave "
      "capillary_coefficient_eV/capillary_constant_eV unset. Exactly one of {Ef, Eb2} or "
      "{capillary_coefficient_eV, capillary_constant_eV} must be given.");
  params.addParam<Real>(
      "Eb2",
      0.0,
      "Dimer binding energy E_b2 [eV]. See Ef for how this derives the capillary law.");
  params.addParam<Real>(
      "capillary_coefficient_eV",
      0.0,
      "Coefficient A [eV] in the capillary formation-energy law E_f(n) = A*n^p + C, given "
      "directly instead of derived from Ef/Eb2 -- the alternative path for a user who already "
      "has a fitted or literature power law (e.g. Soneda's A*n^(2/3) with C = 0, rather than a "
      "monomer/dimer energy pair). The monomer formation energy used for the initial condition "
      "and as E_f(1) in the binding-energy expression becomes A + C. Exactly one of {Ef, Eb2} or "
      "{capillary_coefficient_eV, capillary_constant_eV} must be given.");
  params.addParam<Real>(
      "capillary_constant_eV",
      0.0,
      "Constant C [eV] in the capillary formation-energy law E_f(n) = A*n^p + C. See "
      "capillary_coefficient_eV.");
  params.addParam<Real>(
      "exponent",
      2.0 / 3.0,
      "Exponent p in the capillary formation-energy law E_f(n) = A*n^p + C. Defaults to 2/3 "
      "(surface-area scaling for a 3-D void); set to 0.8 for a typical 2-D loop. Shared by both "
      "the Ef/Eb2 and the direct-coefficient path.");

  params.addParam<MooseEnum>(
      "binding_energy_model",
      MooseEnum("capillary binding_energy_table formation_energy_table", "capillary"),
      "Binding-energy model for this species' emission. 'capillary' (default) always uses the "
      "capillary law derived above. 'binding_energy_table' or 'formation_energy_table' use "
      "direct DFT/MD data for small n (binding_energy_table_eV/formation_energy_table_eV), "
      "falling back to the capillary law for any n beyond the table's range -- the table's own "
      "length is the crossover point, with no separate crossover parameter to set.");
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
      "Z",
      1.0,
      "Dislocation-sink bias Z for this species' monomer (used for the background dislocation "
      "sink, D*Z*dislocation_density). For a loop, Z is ALSO the large-n limit of this species' "
      "own self-absorption bias law below (geometry = loop only).");
  params.addParam<Real>(
      "Z_loop",
      1.0,
      "Coefficient of this species' own self-absorption dislocation bias at the smallest sizes: "
      "Z(n) = Z + [sqrt(b/(8 pi a0)) * Z_loop - Z] * n^(-alpha/2), so Z(1) = sqrt(b/(8 pi a0)) * "
      "Z_loop (about 0.152 * Z_loop for b = a0/sqrt(3)) and Z(n) -> Z as n -> infinity "
      "(geometry = loop only; a void has no bias law).");
  params.addParam<Real>(
      "alpha",
      1.0,
      "Exponent controlling how fast this species' own self-absorption bias approaches Z with "
      "size (geometry = loop only).");

  params.addParam<Real>(
      "Z_loop_opposite",
      1.0,
      "Coefficient of this LOOP's bias for absorbing the OPPOSITE species' monomer at the smallest "
      "sizes, scaled by sqrt(b/(8 pi a0)) as for Z_loop (geometry = loop only, meaningful only "
      "when "
      "both species are present). The large-n limit is the opposite species' own Z; there is "
      "nothing else to supply here.");
  params.addParam<Real>(
      "alpha_opposite",
      1.0,
      "Exponent controlling how fast this loop's opposite-species-absorption bias approaches the "
      "opposite species' Z with size (geometry = loop only).");

  params.addParam<Real>(
      "burgers_vector",
      0.0,
      "Burgers vector b, in m, for a loop's radius and bias laws. Defaults to "
      "lattice_parameter/sqrt(3) (the enclosing `[DefectDynamics]` block's lattice_parameter) when "
      "left at zero (geometry = loop only).");

  params.addParam<MooseEnum>(
      "bias",
      bias,
      "Whether this species' own self-absorption bias (geometry = loop only) uses the full "
      "size-dependent law (Z, Z_loop, alpha) or a single constant value Z at every size (Z only; "
      "Z_loop is then chosen as Z / sqrt(b/(8 pi a0)), see the theory page).");

  params.addParam<unsigned int>(
      "num_cluster_sizes",
      "Maximum physical cluster size (point defects) this species' grouped array represents, "
      "given directly. Void and loop populations routinely need different caps, so this is set "
      "per species, not shared. Exactly one of {num_cluster_sizes, max_cluster_radius_m} must be "
      "given.");
  params.addParam<Real>(
      "max_cluster_radius_m",
      "Maximum physical cluster size, given as a radius [m] instead of a point-defect count; "
      "converted to num_cluster_sizes using this species' own monomer radius law (r_n = "
      "r1*n^(1/3) for a void, r1*n^(1/2) for a loop). Exactly one of {num_cluster_sizes, "
      "max_cluster_radius_m} must be given.");

  // Grouped-array settings of this cluster type. The names are the ones the kernels, the
  // postprocessors and [SoluteClusterDynamics] use. A value left unset falls back to the
  // (deprecated) [Grid] sub-block when there is one, then to these defaults.
  MooseEnum group_spacing("linear geometric", "linear");
  params.addParam<bool>("use_grouping",
                        true,
                        "Whether this cluster type's array stores sizes above explicit_max in "
                        "grouped bins (true) or one unknown per size (false). Without grouping "
                        "every size is an unknown, so it is only practical for small "
                        "num_cluster_sizes.");
  params.addParam<unsigned int>(
      "explicit_max", 100, "Maximum cluster size stored explicitly before the grouped bins.");
  params.addParam<MooseEnum>(
      "group_spacing", group_spacing, "Grouped-bin spacing: linear or geometric.");
  params.addParam<unsigned int>(
      "group_linear_width", 25, "Fixed bin width when group_spacing = linear.");
  params.addParam<Real>(
      "group_geometric_ratio", 1.35, "Bin-width growth factor when group_spacing = geometric.");
  params.addParam<bool>("enforce_group_nonnegative",
                        true,
                        "Whether to limit the grouped L0/L1 reconstruction so bin endpoint "
                        "concentrations stay nonnegative (up to the tolerance factor).");
  params.addParam<Real>("group_nonnegative_tolerance_factor",
                        0.0,
                        "Relaxation factor for the nonnegativity limit: endpoint concentrations "
                        "may go down to -factor * L0 before being limited.");

  params.addParam<Real>(
      "resolvable_radius",
      0.0,
      "If positive, adds a number-density and average-radius postprocessor pair "
      "(<name>_resolvable_density, <name>_resolvable_avg_radius) for this "
      "species, counting only clusters at or above this radius [m] -- the smallest cluster size "
      "resolvable by whatever characterization technique the input is being compared against "
      "(e.g. TEM). Left at zero (the default), no such postprocessors are added for this "
      "species.");

  params.addParam<bool>(
      "output_size_distribution",
      false,
      "If true, adds a ClusterSizeDistributionGrouped vector postprocessor for this species, "
      "reporting its full per-size concentration (and radius) distribution at the enclosing "
      "DefectDynamics block's size_distribution_times. Left false (the default), no such output "
      "is added for this species.");

  params.addParam<Real>(
      "cascade_f2",
      0.0,
      "Fraction of this species' generation rate seeded directly as dimers (cluster size 2) "
      "instead of monomers (Pokor et al. 2004, Eq 1). A nonzero cascade_f2, cascade_f3 or "
      "cascade_f4 turns cascade seeding on for this species, and the monomer then receives the "
      "rest of the generation rate. Cascade seeding needs sizes 2 to 4 to be explicit array "
      "components, so the species must have explicit_max >= 4. See cascade_f3, cascade_f4.");
  params.addParam<Real>(
      "cascade_f3",
      0.0,
      "Fraction of this species' generation rate seeded directly as cluster size 3. See "
      "cascade_f2.");
  params.addParam<Real>(
      "cascade_f4",
      0.0,
      "Fraction of this species' generation rate seeded directly as cluster size 4. See "
      "cascade_f2.");

  return params;
}

namespace
{
/// Largest num_cluster_sizes for a species that does not group its array.
constexpr unsigned int max_ungrouped_cluster_sizes = 3000;

/// Largest cluster size a radius may be converted to. The kernels build arrays of this length on
/// every call, so a size near the limit of an unsigned int would be out of reach anyway.
constexpr Real max_cluster_size_from_radius = 1.0e9;

bool
checkCapillaryInputPair(const InputParameters & params)
{
  const bool has_Ef_Eb2 =
      params.isParamSetByUser("Ef") || params.isParamSetByUser("Eb2");
  const bool has_direct = params.isParamSetByUser("capillary_coefficient_eV") ||
                         params.isParamSetByUser("capillary_constant_eV");
  if (has_Ef_Eb2 && has_direct)
    mooseError("DefectSpeciesAction cannot set both {Ef, Eb2} and {capillary_coefficient_eV, "
              "capillary_constant_eV}; use exactly one of the two ways to specify the capillary "
              "formation-energy law.");
  if (!has_Ef_Eb2 && !has_direct)
    mooseError("DefectSpeciesAction requires either {Ef, Eb2} or {capillary_coefficient_eV, "
              "capillary_constant_eV} to specify the capillary formation-energy law.");
  return has_direct;
}

bool
checkClusterSizeInputPair(const InputParameters & params)
{
  const bool has_count = params.isParamSetByUser("num_cluster_sizes");
  const bool has_radius = params.isParamSetByUser("max_cluster_radius_m");
  if (has_count && has_radius)
    mooseError("DefectSpeciesAction cannot set both num_cluster_sizes and "
              "max_cluster_radius_m; use exactly one of the two ways to specify this species' "
              "maximum cluster size.");
  if (!has_count && !has_radius)
    mooseError("DefectSpeciesAction requires either num_cluster_sizes or max_cluster_radius_m "
              "to specify this species' maximum cluster size.");
  return has_radius;
}
} // namespace

DefectSpeciesAction::DefectSpeciesAction(const InputParameters & params)
  : Action(params),
    _uses_direct_capillary_coefficients(checkCapillaryInputPair(params)),
    _uses_direct_cluster_size_count(!checkClusterSizeInputPair(params))
{
  const auto geometry = getParam<MooseEnum>("geometry");
  if (geometry == "void" && parameters().isParamSetByUser("Z_loop"))
    paramWarning("Z_loop", "Z_loop has no effect for geometry = void: a void has no "
                           "dislocation-bias absorption law.");
  if (geometry == "void" && parameters().isParamSetByUser("alpha"))
    paramWarning(
        "alpha", "alpha has no effect for geometry = void: a void has no dislocation-bias "
                "absorption law.");
  if (geometry == "void" && parameters().isParamSetByUser("burgers_vector"))
    paramWarning("burgers_vector",
                "burgers_vector has no effect for geometry = void: only a loop's radius and "
                "bias laws use it.");

  if (parameters().isParamSetByUser("max_cluster_radius_m") &&
      getParam<Real>("max_cluster_radius_m") <= 0.0)
    paramError("max_cluster_radius_m", "must be positive.");
  if (getParam<Real>("resolvable_radius") < 0.0)
    paramError("resolvable_radius", "must not be negative; leave it at zero to turn it off.");

  const Real cascade_f2 = getParam<Real>("cascade_f2");
  const Real cascade_f3 = getParam<Real>("cascade_f3");
  const Real cascade_f4 = getParam<Real>("cascade_f4");
  if (cascade_f2 < 0.0 || cascade_f3 < 0.0 || cascade_f4 < 0.0)
    mooseError("DefectSpeciesAction requires cascade_f2, cascade_f3, and cascade_f4 to be "
              "nonnegative.");
  if (cascade_f2 + cascade_f3 + cascade_f4 > 1.0)
    mooseError("DefectSpeciesAction requires cascade_f2 + cascade_f3 + cascade_f4 <= 1, so the "
              "monomer's own remaining fraction of the generation rate does not go negative.");
}

DefectSpeciesAction::Grouping
DefectSpeciesAction::grouping(const DefectDynamicsGridAction * grid) const
{
  // A value this species set itself wins; otherwise the [Grid] sub-block's, if there is one;
  // otherwise this species' own default (identical to the [Grid] defaults).
  auto own = [this](const std::string & name) { return parameters().isParamSetByUser(name); };

  Grouping g;
  g.use_grouping = getParam<bool>("use_grouping");
  g.explicit_max = own("explicit_max") || !grid ? getParam<unsigned int>("explicit_max")
                                                : grid->getParam<unsigned int>("explicit_max");
  g.spacing = own("group_spacing") || !grid ? getParam<MooseEnum>("group_spacing")
                                            : grid->getParam<MooseEnum>("spacing");
  g.linear_width = own("group_linear_width") || !grid
                       ? getParam<unsigned int>("group_linear_width")
                       : grid->getParam<unsigned int>("linear_width");
  g.geometric_ratio = own("group_geometric_ratio") || !grid
                          ? getParam<Real>("group_geometric_ratio")
                          : grid->getParam<Real>("geometric_ratio");
  g.enforce_nonnegative = own("enforce_group_nonnegative") || !grid
                              ? getParam<bool>("enforce_group_nonnegative")
                              : grid->getParam<bool>("enforce_nonnegative");
  g.tolerance_factor = own("group_nonnegative_tolerance_factor") || !grid
                           ? getParam<Real>("group_nonnegative_tolerance_factor")
                           : grid->getParam<Real>("tolerance_factor");
  return g;
}

ClusterGroupingLayout
DefectSpeciesAction::layout(unsigned int num_cluster_sizes,
                            const DefectDynamicsGridAction * grid) const
{
  const auto g = grouping(grid);
  return ClusterGroupingLayout(num_cluster_sizes,
                               g.explicit_max,
                               g.use_grouping,
                               g.spacing.getEnum<ClusterGroupingSpacing>(),
                               g.linear_width,
                               g.geometric_ratio);
}

void
DefectSpeciesAction::monomerRadiusLaw(
    Real atomic_volume, Real lattice_parameter, Real & r1, Real & radius_exponent) const
{
  if (getParam<MooseEnum>("geometry") == "void")
  {
    r1 = std::cbrt(3.0 * atomic_volume / (4.0 * libMesh::pi));
    radius_exponent = 1.0 / 3.0;
  }
  else
  {
    const Real burgers_vector = getParam<Real>("burgers_vector") != 0.0
                                     ? getParam<Real>("burgers_vector")
                                     : lattice_parameter / std::sqrt(3.0);
    r1 = std::sqrt(atomic_volume / (libMesh::pi * burgers_vector));
    radius_exponent = 0.5;
  }
}

unsigned int
DefectSpeciesAction::clusterSizeFromRadius(Real radius, Real r1, Real radius_exponent)
{
  // Converting a double that does not fit an unsigned int is undefined behavior (and, for a large
  // radius, ends in an allocation of billions of sizes), so reject it with a message instead.
  const Real size = std::ceil(std::pow(radius / r1, 1.0 / radius_exponent));
  if (!std::isfinite(size) || size < 1.0 || size > max_cluster_size_from_radius)
    ::mooseError("A radius of ",
                 radius,
                 " m corresponds to ",
                 size,
                 " point defects with a monomer radius of ",
                 r1,
                 " m; it must correspond to between 1 and ",
                 max_cluster_size_from_radius,
                 ".");
  return static_cast<unsigned int>(size);
}

unsigned int
DefectSpeciesAction::resolvedNumClusterSizes(Real atomic_volume, Real lattice_parameter) const
{
  unsigned int n;
  if (_uses_direct_cluster_size_count)
    n = getParam<unsigned int>("num_cluster_sizes");
  else
  {
    Real r1, radius_exponent;
    monomerRadiusLaw(atomic_volume, lattice_parameter, r1, radius_exponent);
    n = clusterSizeFromRadius(getParam<Real>("max_cluster_radius_m"), r1, radius_exponent);
  }

  // The monomer-row balance reads the dimer, so a species needs at least sizes 1 and 2.
  if (n < 2)
    paramError(_uses_direct_cluster_size_count ? "num_cluster_sizes" : "max_cluster_radius_m",
               "a species needs at least 2 cluster sizes, but this one has ",
               n,
               ".");

  // Without grouping every size is an unknown, and both the direct solve and the assembly storage
  // grow much faster than linearly with the unknown count, so a very large ungrouped species would
  // run for hours or exhaust memory without any message.
  if (!getParam<bool>("use_grouping") && n > max_ungrouped_cluster_sizes)
    paramError("use_grouping",
               "use_grouping = false gives every cluster size its own unknown, so a species may "
               "have at most ",
               max_ungrouped_cluster_sizes,
               " sizes (this one has ",
               n,
               "). Set use_grouping = true, or lower the maximum cluster size.");
  return n;
}
