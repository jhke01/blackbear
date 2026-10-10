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

/**
 * The [DefectDynamics/Grid] sub-block: the grouped-tail SHAPE (explicit_max, spacing,
 * linear_width, geometric_ratio, enforce_nonnegative, tolerance_factor) shared by every species
 * added under the enclosing [DefectDynamics] block. The maximum cluster SIZE itself
 * (num_cluster_sizes) is deliberately not described here -- each species sets its own, directly
 * or via a maximum cluster radius, on its own DefectSpeciesAction block, since void and loop
 * populations routinely need very different size caps. Holds no logic of its own; it exists only
 * so DefectDynamicsAction can look it up (ActionWarehouse::getActions<DefectDynamicsGridAction>())
 * and read its validated parameters while building the final variables, ICs, kernels, and
 * postprocessors.
 */
class DefectDynamicsGridAction : public Action
{
public:
  static InputParameters validParams();

  DefectDynamicsGridAction(const InputParameters & params);

  virtual void act() override {}

  /// The grouping layout described by this block's shape parameters, for the given (per-species)
  /// maximum cluster size.
  ClusterGroupingLayout layout(unsigned int num_cluster_sizes) const;
};
