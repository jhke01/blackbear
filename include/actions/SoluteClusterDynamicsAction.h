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

class DefectSpeciesAction;

/**
 * [SoluteClusterDynamics]: a user-facing front end to the solute-clustering ("solute cluster
 * dynamics") kernels, GenericClusterDynamicsGroupedNodalKernel and
 * GenericClusterDynamicsNodalKernel. One block builds the solute species' array variable, initial
 * condition, time-derivative/ dynamics kernel, and monomer/total-density postprocessors, as a thin
 * pass-through of the kernel's own physical parameters -- unlike [DefectDynamics], the solute
 * kernel already derives its own rate coefficients internally, so this action has little physics of
 * its own to compute.
 *
 * This block's main value-add is automatic point-defect coupling: when a [DefectDynamics] block
 * is ALSO present in the same input, this action looks for its vacancy species (identified by
 * geometry = void, falling back to a species literally named "vacancy") and, when found, wires
 * vacancy_coupled_var/vacancy_formation_energy_eV onto the solute kernel automatically. No
 * separate "activate coupling" switch exists; declaring both blocks together is how a user
 * studying solute clusters opts into concurrent point-defect clustering and its effect on solute
 * diffusivity. A [SoluteClusterDynamics] block with no [DefectDynamics] block present builds the
 * solute species exactly as it would otherwise, with the coupling simply inactive.
 */
class SoluteClusterDynamicsAction : public Action
{
public:
  static InputParameters validParams();

  SoluteClusterDynamicsAction(const InputParameters & params);

  virtual void act() override;

protected:
  void actAddVariables();
  void actAddFunctions();
  void actAddICs();
  void actAddNodalKernels();
  void actAddPostprocessors();

  /// The array variable name this block builds, "<block name>_clusters".
  std::string variableName() const;
  /// The initial-condition function name this block builds, "<block name>_index_fn".
  std::string functionName() const;

  /// Looks for a [DefectDynamics] block's vacancy species to couple to: nullptr if no
  /// [DefectDynamics] block is present (the normal, uncoupled case). Errors if more than one
  /// [DefectDynamics] block exists, or if a single block's candidate vacancy species (by
  /// geometry = void, falling back to a species named "vacancy") is ambiguous.
  const DefectSpeciesAction * findVacancySpecies() const;
};
