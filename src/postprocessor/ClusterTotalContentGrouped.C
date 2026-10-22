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

#include "ClusterTotalContentGrouped.h"

#include <limits>

registerMooseObject("BlackBearApp", ClusterTotalContentGrouped);

InputParameters
ClusterTotalContentGrouped::validParams()
{
  InputParameters params = NodalPostprocessor::validParams();
  params.set<bool>("unique_node_execute") = true;
  MooseEnum group_spacing("linear geometric", "linear");
  params.addClassDescription(
      "Computes the total point-defect content (sum of n * C_n) represented by an explicit and "
      "grouped L0/L1 cluster dynamics array variable.");
  params.addRequiredCoupledVar("clusters", "Grouped cluster array concentration variable.");
  params.addRequiredParam<unsigned int>(
      "num_cluster_sizes", "Maximum physical cluster size represented by the grouped array.");
  params.addParam<bool>(
      "use_grouping", true, "Whether the coupled array uses grouped-tail storage.");
  params.addParam<unsigned int>(
      "explicit_max", 100, "Maximum physical cluster size stored explicitly before grouped bins.");
  params.addParam<MooseEnum>(
      "group_spacing", group_spacing, "Grouped-tail spacing mode: linear or geometric.");
  params.addParam<unsigned int>(
      "group_linear_width", 25, "Fixed grouped-bin width when group_spacing = linear.");
  params.addParam<Real>("group_geometric_ratio",
                        1.35,
                        "Geometric grouped-bin growth factor when group_spacing = geometric.");
  params.addParam<bool>(
      "enforce_group_nonnegative",
      true,
      "Whether to limit grouped L0/L1 reconstruction to keep bin endpoint concentrations "
      "nonnegative (up to the tolerance factor).");
  params.addParam<Real>(
      "group_nonnegative_tolerance_factor",
      0.0,
      "Relaxation factor for grouped nonnegativity enforcement. Endpoint concentrations are "
      "allowed down to -factor * L0 before limiting.");
  params.addRangeCheckedParam<unsigned int>(
      "n_minimum",
      1,
      "n_minimum > 0",
      "Minimum physical cluster size to include in the total content (default: 1, includes "
      "monomers).");
  params.addParam<Real>("atomic_volume", 1.0, "Atomic volume [m^3] for output conversion.");
  return params;
}

ClusterTotalContentGrouped::ClusterTotalContentGrouped(const InputParameters & parameters)
  : NodalPostprocessor(parameters),
    _clusters(coupledArrayValue("clusters")),
    _n_minimum(getParam<unsigned int>("n_minimum")),
    _atomic_volume(getParam<Real>("atomic_volume")),
    _layout(getParam<unsigned int>("num_cluster_sizes"),
            getParam<unsigned int>("explicit_max"),
            getParam<bool>("use_grouping"),
            getParam<MooseEnum>("group_spacing").getEnum<ClusterGroupingSpacing>(),
            getParam<unsigned int>("group_linear_width"),
            getParam<Real>("group_geometric_ratio")),
    _enforce_group_nonnegative(getParam<bool>("enforce_group_nonnegative")),
    _group_nonnegative_tolerance_factor(getParam<Real>("group_nonnegative_tolerance_factor")),
    _sum(0.0),
    _count(0.0)
{
  if (_atomic_volume <= 0.0)
    mooseError("ClusterTotalContentGrouped requires atomic_volume > 0.");
  if (_group_nonnegative_tolerance_factor < 0.0)
    mooseError("ClusterTotalContentGrouped requires group_nonnegative_tolerance_factor >= 0.");
  if (getArrayVar("clusters", 0)->count() != _layout.componentCount())
    paramError("clusters",
               "The grouping layout requires ",
               _layout.componentCount(),
               " array components, but the coupled variable has ",
               getArrayVar("clusters", 0)->count(),
               ". The grouping parameters must match those of the cluster variable.");
}

Real
ClusterTotalContentGrouped::concentrationAt(unsigned int n) const
{
  return _layout.concentration(
      _clusters[_qp], n, _enforce_group_nonnegative, _group_nonnegative_tolerance_factor);
}

bool
ClusterTotalContentGrouped::groupedBinNonnegative(const ClusterGroupingBin & bin) const
{
  if (!_enforce_group_nonnegative)
    return false;

  const Real l0 = _clusters[_qp](bin.l0_component);
  if (l0 <= 0.0)
    return false;

  const Real l1 = _clusters[_qp](bin.l1_component);
  const Real left_dx = static_cast<Real>(bin.start) - bin.mean_x;
  const Real right_dx = static_cast<Real>(bin.end) - bin.mean_x;
  Real lower = -std::numeric_limits<Real>::infinity();
  Real upper = std::numeric_limits<Real>::infinity();
  if (right_dx > 0.0)
    lower = -l0 / right_dx;
  if (left_dx < 0.0)
    upper = l0 / (-left_dx);

  return l1 >= lower && l1 <= upper;
}

void
ClusterTotalContentGrouped::initialize()
{
  _sum = 0.0;
  _count = 0.0;
}

void
ClusterTotalContentGrouped::execute()
{
  Real node_total = 0.0;
  const unsigned int explicit_upper = std::min(_layout.explicitMax(), _layout.numClusterSizes());
  for (unsigned int n = _n_minimum; n <= explicit_upper; ++n)
    node_total += n * concentrationAt(n);

  for (const auto & bin : _layout.groups())
  {
    if (bin.end < _n_minimum)
      continue;

    // With enforce_group_nonnegative and L0 <= 0 the layout reconstructs every size in the bin as
    // exactly zero (reconstructionCoefficients() returns a0 = a1 = 0), so the bin adds only +-0 to
    // the sums. Skipping it avoids a per-size loop over the empty tail on every step.
    if (_enforce_group_nonnegative && _clusters[_qp](bin.l0_component) <= 0.0)
      continue;

    // For a complete, unclipped bin, sum_n n*C_n = L0*width*mean_x + L1*width*sigma2: the L1 term
    // reduces to width*sigma2 because sum_n (n-mean_x)*n = sum_n (n-mean_x)^2 (the mean_x*sum_n
    // (n-mean_x) cross term vanishes by definition of mean_x).
    if (bin.start >= _n_minimum && groupedBinNonnegative(bin))
    {
      const Real l0 = _clusters[_qp](bin.l0_component);
      const Real l1 = _clusters[_qp](bin.l1_component);
      node_total += bin.width * (l0 * bin.mean_x + l1 * bin.sigma2);
      continue;
    }

    for (unsigned int n = std::max(_n_minimum, bin.start); n <= bin.end; ++n)
      node_total += n * concentrationAt(n);
  }

  _sum += node_total;
  _count += 1.0;
}

void
ClusterTotalContentGrouped::finalize()
{
  gatherSum(_sum);
  gatherSum(_count);
}

void
ClusterTotalContentGrouped::threadJoin(const UserObject & y)
{
  const auto & rhs = static_cast<const ClusterTotalContentGrouped &>(y);
  _sum += rhs._sum;
  _count += rhs._count;
}

Real
ClusterTotalContentGrouped::getValue() const
{
  return _count > 0.0 ? _sum / (_count * _atomic_volume) : 0.0;
}
