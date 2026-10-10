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
 * The [DefectDynamics/Grid] sub-block: the single shared grouped-tail layout used by every
 * species added under the enclosing [DefectDynamics] block. Holds no logic of its own; it exists
 * only so DefectDynamicsAction can look it up (ActionWarehouse::getActions<DefectDynamicsGridAction>())
 * and read its validated parameters while building the final variables, ICs, kernels, and
 * postprocessors.
 */
class DefectDynamicsGridAction : public Action
{
public:
  static InputParameters validParams();

  DefectDynamicsGridAction(const InputParameters & params);

  virtual void act() override {}

  /// The grouping layout described by this block's parameters.
  ClusterGroupingLayout layout() const;
};
