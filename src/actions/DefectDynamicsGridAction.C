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
      "The grouped-tail layout shared by every species in the enclosing DefectDynamics block.");
  params.addRequiredParam<unsigned int>(
      "num_cluster_sizes", "Maximum physical cluster size represented by the grouped array.");
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
DefectDynamicsGridAction::layout() const
{
  return ClusterGroupingLayout(getParam<unsigned int>("num_cluster_sizes"),
                               getParam<unsigned int>("explicit_max"),
                               true,
                               getParam<MooseEnum>("spacing").getEnum<ClusterGroupingSpacing>(),
                               getParam<unsigned int>("linear_width"),
                               getParam<Real>("geometric_ratio"));
}
