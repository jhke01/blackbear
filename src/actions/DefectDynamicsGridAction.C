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

#include "DefectDynamicsGridAction.h"

registerMooseAction("BlackBearApp", DefectDynamicsGridAction, "meta_action");

InputParameters
DefectDynamicsGridAction::validParams()
{
  InputParameters params = Action::validParams();
  MooseEnum group_spacing("linear geometric", "linear");
  params.addClassDescription(
      "The grouped-tail shape shared by every species in the enclosing DefectDynamics block. "
      "The maximum cluster size itself is set per species (DefectSpeciesAction's "
      "num_cluster_sizes or max_cluster_radius_m), not here.");
  params.addParam<unsigned int>(
      "explicit_max", 100, "Maximum physical cluster size stored explicitly before grouped bins.");
  params.addParam<MooseEnum>(
      "spacing", group_spacing, "Grouped-tail spacing mode: linear or geometric.");
  params.addParam<unsigned int>(
      "linear_width", 25, "Fixed grouped-bin width when spacing = linear.");
  params.addParam<Real>(
      "geometric_ratio", 1.35, "Geometric grouped-bin growth factor when spacing = geometric.");
  params.addParam<bool>(
      "enforce_nonnegative",
      true,
      "Whether to limit grouped L0/L1 reconstruction to keep bin endpoint concentrations "
      "nonnegative (up to the tolerance factor).");
  params.addParam<Real>(
      "tolerance_factor",
      0.0,
      "Relaxation factor for grouped nonnegativity enforcement. Endpoint concentrations are "
      "allowed down to -factor * L0 before limiting.");
  return params;
}

DefectDynamicsGridAction::DefectDynamicsGridAction(const InputParameters & params) : Action(params)
{
}

ClusterGroupingLayout
DefectDynamicsGridAction::layout(unsigned int num_cluster_sizes) const
{
  return ClusterGroupingLayout(num_cluster_sizes,
                               getParam<unsigned int>("explicit_max"),
                               true,
                               getParam<MooseEnum>("spacing").getEnum<ClusterGroupingSpacing>(),
                               getParam<unsigned int>("linear_width"),
                               getParam<Real>("geometric_ratio"));
}
