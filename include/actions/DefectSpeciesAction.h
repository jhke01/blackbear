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

#pragma once

#include "Action.h"
#include "MooseEnum.h"

/**
 * One defect species sub-block under [DefectDynamics] (e.g. [vacancy], [interstitial]). The
 * sub-block's own name is just a label; "geometry" selects void (3-D, finite-sink-strength
 * absorption) or loop (2-D, dislocation-bias absorption) behavior, matching
 * GenericVacancyVoidGroupedNodalKernel/GenericInterstitialLoopGroupedNodalKernel respectively.
 *
 * This action holds no logic of its own -- act() is a no-op. It exists purely so the enclosing
 * DefectDynamicsAction can find every species block via ActionWarehouse::getActions<DefectSpeciesAction>()
 * and read its validated parameters when building the final variables, ICs, kernels, and
 * postprocessors. One or two species sub-blocks are expected; with two, each species'
 * "opposite_*" kernel parameters (bias, diffusivity, grid) are derived entirely from the other
 * species' own block, so the user never types anything named "opposite".
 */
class DefectSpeciesAction : public Action
{
public:
  static InputParameters validParams();

  DefectSpeciesAction(const InputParameters & params);

  virtual void act() override {}

  /// Whether this species specified Ef/Eb2 (the default path) or capillary_coefficient_eV/
  /// capillary_constant_eV directly (the alternative path). Exactly one of the two pairs must be
  /// set; validated in the constructor.
  bool usesDirectCapillaryCoefficients() const { return _uses_direct_capillary_coefficients; }

  /// Whether this species specified num_cluster_sizes directly, or max_cluster_radius_m instead.
  /// Exactly one of the two must be set; validated in the constructor.
  bool usesDirectClusterSizeCount() const { return _uses_direct_cluster_size_count; }

  /// This species' own maximum cluster size, resolved from whichever of num_cluster_sizes or
  /// max_cluster_radius_m was given -- the latter converted via this species' own monomer radius
  /// law (void or loop, using the enclosing [DefectDynamics] block's atomic_volume and
  /// lattice_parameter).
  unsigned int resolvedNumClusterSizes(Real atomic_volume, Real lattice_parameter) const;

  /// r1 (the monomer radius) and the power-law exponent in r_n = r1*n^exponent for this species'
  /// own geometry (void: exponent 1/3; loop: exponent 1/2, using this species' own burgers_vector,
  /// defaulting to lattice_parameter/sqrt(3)). Shared by resolvedNumClusterSizes() and by
  /// DefectDynamicsAction's resolvable_radius -> n_minimum conversion, so both use the identical
  /// radius law for this species rather than two separately maintained copies of it.
  void monomerRadiusLaw(Real atomic_volume, Real lattice_parameter, Real & r1, Real & radius_exponent) const;

  /// Inverts r_n = r1*n^radius_exponent to the smallest integer cluster size n whose radius is at
  /// least `radius`. Shared by resolvedNumClusterSizes() (max_cluster_radius_m) and
  /// DefectDynamicsAction's resolvable_radius (n_minimum) -- the same inversion, two different
  /// physical radii.
  static unsigned int clusterSizeFromRadius(Real radius, Real r1, Real radius_exponent);

protected:
  const bool _uses_direct_capillary_coefficients;
  const bool _uses_direct_cluster_size_count;
};
