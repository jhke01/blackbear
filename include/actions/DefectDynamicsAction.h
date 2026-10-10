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
#include "ClusterGroupingLayout.h"

class DefectDynamicsGridAction;
class DefectRecoveryAction;
class DefectSpeciesAction;

/**
 * [DefectDynamics]: a user-facing, species-organized front end to the grouped defect-dynamics
 * kernels (GenericVacancyVoidGroupedNodalKernel, GenericInterstitialLoopGroupedNodalKernel,
 * PointDefectRecombinationNodalKernel). The block declares the physics shared by every species
 * (temperature, atomic_volume, lattice_parameter, cascade generation, dislocation sink,
 * recombination), one [Grid] sub-block (DefectDynamicsGridAction) describing the grouped-tail
 * SHAPE shared by every species, and one or two named species sub-blocks (DefectSpeciesAction,
 * e.g. [vacancy], [interstitial]) describing each species' own diffusivity, capillary energy law,
 * (for a loop) dislocation bias, and its own maximum cluster size (void and loop populations
 * routinely need very different size caps, so this is not shared via [Grid]).
 *
 * This action builds every final object -- the grouped array variables, their initial
 * conditions, the self/opposite-absorption kernels, the recombination kernels (when two species
 * are present), and the monomer/resolvable-radius postprocessors -- deriving every quantity the
 * low-level kernels call "opposite_*" from the OTHER species' own block, and deriving the
 * dislocation sink, recombination rate, capillary coefficients, and monomer initial condition
 * from the physics parameters directly, rather than asking the user to type any of them.
 *
 * Per BlackBear's existing convention (ClusterGroupingLayout, GroupedClusterVariableAction), the
 * low-level kernels and the [GroupedVariables]/[NodalKernels] input idiom remain available
 * unchanged for a user who needs an asymmetric grid, a single species, or direct control over a
 * parameter this action does not expose.
 */
class DefectDynamicsAction : public Action
{
public:
  static InputParameters validParams();

  DefectDynamicsAction(const InputParameters & params);

  virtual void act() override;

protected:
  /// This block's [Grid] sub-block, or nullptr if none was given (an error, checked in act()).
  const DefectDynamicsGridAction * grid() const;
  /// This block's optional [Recovery] sub-block, or nullptr if none was given (zero is the
  /// normal, default-off case; more than one is an error).
  const DefectRecoveryAction * recovery() const;
  /// This block's species sub-blocks (one or two), sorted by block name.
  std::vector<const DefectSpeciesAction *> species() const;

  /// 4*exp(-Ef/kT), the detailed-balance monomer concentration for a species with formation
  /// energy Ef [eV] at this block's temperature, using the SAME kB/eV constants as the kernels'
  /// own anonymous-namespace values so the derived number is bit-identical to a user-computed one
  /// using those same constants, not merely numerically close.
  Real monomerConcentration(Real Ef_eV) const;

  /// A*n^exponent + C, the capillary formation-energy coefficient/constant pair derived from
  /// (Ef, Eb2, exponent): A = (Ef - Eb2)/(2^exponent - 1), C = Ef - A.
  static void capillaryCoefficients(Real Ef_eV, Real Eb2_eV, Real exponent, Real & A, Real & C);

  /// Resolves a species' capillary A, C, and monomer formation energy (E_f(1) = A + C) from
  /// whichever of the two mutually exclusive input pairs it set: {Ef, Eb2} (derived via
  /// capillaryCoefficients) or {capillary_coefficient_eV, capillary_constant_eV} (used directly).
  static void resolveCapillaryCoefficients(
      const DefectSpeciesAction & sp, Real & A, Real & C, Real & formation_energy_eV);

  void actAddVariables();
  void actAddFunctions();
  void actAddICs();
  void actAddNodalKernels();
  void actAddPostprocessors();

  /// Builds the array variable name for a species block, e.g. "vacancy_clusters" for a block
  /// named "vacancy". Shared across actAddVariables/actAddFunctions/actAddICs/actAddNodalKernels/
  /// actAddPostprocessors so every object referencing a species agrees on its variable's name.
  static std::string variableName(const DefectSpeciesAction & sp);
  static std::string functionName(const DefectSpeciesAction & sp);
};
