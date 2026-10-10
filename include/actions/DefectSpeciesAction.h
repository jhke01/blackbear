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

protected:
  const bool _uses_direct_capillary_coefficients;
};
