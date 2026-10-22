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

#include "ClusterSizeDistributionGrouped.h"

#include <cmath>

registerMooseObject("BlackBearApp", ClusterSizeDistributionGrouped);

InputParameters
ClusterSizeDistributionGrouped::validParams()
{
  InputParameters params = NodalVectorPostprocessor::validParams();
  params.set<bool>("unique_node_execute") = true;
  MooseEnum group_spacing("linear geometric", "linear");
  params.addClassDescription(
      "Reports the full per-size cluster concentration distribution (one row per physical "
      "cluster size n) from a grouped L0/L1 cluster dynamics array variable, reconstructing "
      "grouped (large-n) sizes the same way the grouped kernels do internally.");
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
  params.addParam<Real>(
      "r1", 1.0, "Monomer radius [m] in the radius law r_n = r1*n^radius_exponent.");
  params.addParam<Real>(
      "radius_exponent",
      1.0 / 3.0,
      "Exponent in the radius law r_n = r1*n^radius_exponent (1/3 for a void, 1/2 for a loop).");
  return params;
}

ClusterSizeDistributionGrouped::ClusterSizeDistributionGrouped(const InputParameters & parameters)
  : NodalVectorPostprocessor(parameters),
    _clusters(coupledArrayValue("clusters")),
    _r1(getParam<Real>("r1")),
    _radius_exponent(getParam<Real>("radius_exponent")),
    _layout(getParam<unsigned int>("num_cluster_sizes"),
            getParam<unsigned int>("explicit_max"),
            getParam<bool>("use_grouping"),
            getParam<MooseEnum>("group_spacing").getEnum<ClusterGroupingSpacing>(),
            getParam<unsigned int>("group_linear_width"),
            getParam<Real>("group_geometric_ratio")),
    _enforce_group_nonnegative(getParam<bool>("enforce_group_nonnegative")),
    _group_nonnegative_tolerance_factor(getParam<Real>("group_nonnegative_tolerance_factor")),
    _n(declareVector("n")),
    _concentration(declareVector("concentration")),
    _radius(declareVector("radius"))
{
  const unsigned int num_cluster_sizes = _layout.numClusterSizes();
  _n_table.resize(num_cluster_sizes);
  _radius_table.resize(num_cluster_sizes);
  for (unsigned int n = 1; n <= num_cluster_sizes; ++n)
  {
    _n_table[n - 1] = n;
    _radius_table[n - 1] = _r1 * std::pow(static_cast<Real>(n), _radius_exponent);
  }

  if (_group_nonnegative_tolerance_factor < 0.0)
    mooseError(
        "ClusterSizeDistributionGrouped requires group_nonnegative_tolerance_factor >= 0.");
  if (getArrayVar("clusters", 0)->count() != _layout.componentCount())
    paramError("clusters",
               "The grouping layout requires ",
               _layout.componentCount(),
               " array components, but the coupled variable has ",
               getArrayVar("clusters", 0)->count(),
               ". The grouping parameters must match those of the cluster variable.");
}

void
ClusterSizeDistributionGrouped::initialize()
{
  const unsigned int num_cluster_sizes = _layout.numClusterSizes();
  _n.resize(num_cluster_sizes);
  _concentration.resize(num_cluster_sizes);
  _radius.resize(num_cluster_sizes);
}

void
ClusterSizeDistributionGrouped::execute()
{
  // This model is always run on a single-node (0D) mesh in every current production and test
  // input, so one node's own execute() call covers the entire distribution; there is no
  // cross-node merge to perform in finalize() (NodalVectorPostprocessor::finalize() defaults to a
  // no-op, which is correct here).
  for (unsigned int n = 1; n <= _layout.numClusterSizes(); ++n)
  {
    const unsigned int i = n - 1;
    _concentration[i] = _layout.concentration(
        _clusters[_qp], n, _enforce_group_nonnegative, _group_nonnegative_tolerance_factor);
  }
  _n = _n_table;
  _radius = _radius_table;
}

void
ClusterSizeDistributionGrouped::threadJoin(const UserObject & /*y*/)
{
  // No-op: this model is always run on a single-node (0D) mesh (see execute()), so there is
  // exactly one node's worth of data and nothing to merge across threads.
}
