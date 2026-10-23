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

#include "ClusterContinuumDistribution.h"

#include <cmath>

registerMooseObject("BlackBearApp", ClusterContinuumDistribution);

InputParameters
ClusterContinuumDistribution::validParams()
{
  InputParameters params = NodalVectorPostprocessor::validParams();
  params.set<bool>("unique_node_execute") = true;
  params.addClassDescription(
      "Reports the size distribution of a cluster variable that is discrete up to switch_size and "
      "a Fokker-Planck continuum above it: one row per discrete size and per continuum cell.");
  params.addRequiredCoupledVar("clusters", "Cluster array concentration variable.");
  params.addRequiredParam<unsigned int>("num_cluster_sizes",
                                        "Largest cluster size represented, N.");
  params.addRequiredParam<unsigned int>(
      "switch_size", "Largest cluster size that is an individual unknown; larger sizes are cells.");
  params.addParam<Real>("continuum_initial_width", 1.0, "Width of the first continuum cell.");
  params.addParam<Real>(
      "continuum_width_ratio", 1.04, "Growth factor of the cell width from one cell to the next.");
  params.addParam<Real>("r1", 1.0, "Monomer radius [m] in r_n = r1*n^radius_exponent.");
  params.addParam<Real>("radius_exponent",
                        0.5,
                        "Exponent in r_n = r1*n^radius_exponent: 1/2 for loops, 1/3 for voids.");
  return params;
}

ClusterContinuumDistribution::ClusterContinuumDistribution(const InputParameters & parameters)
  : NodalVectorPostprocessor(parameters),
    _clusters(coupledArrayValue("clusters")),
    _r1(getParam<Real>("r1")),
    _radius_exponent(getParam<Real>("radius_exponent")),
    _layout(getParam<unsigned int>("num_cluster_sizes"),
            getParam<unsigned int>("switch_size"),
            getParam<Real>("continuum_initial_width"),
            getParam<Real>("continuum_width_ratio")),
    _n(declareVector("n")),
    _width(declareVector("width")),
    _concentration(declareVector("concentration")),
    _radius(declareVector("radius"))
{
  if (getArrayVar("clusters", 0)->count() != _layout.componentCount())
    paramError("clusters",
               "The continuum layout requires ",
               _layout.componentCount(),
               " array components, but the coupled variable has ",
               getArrayVar("clusters", 0)->count(),
               ". The layout parameters must match those of the cluster variable.");
}

void
ClusterContinuumDistribution::initialize()
{
  const unsigned int rows = _layout.switchSize() + _layout.numCells();
  _n.resize(rows);
  _width.resize(rows);
  _concentration.resize(rows);
  _radius.resize(rows);
}

void
ClusterContinuumDistribution::execute()
{
  const RealEigenVector & state = _clusters[_qp];
  auto radius = [this](Real size) { return _r1 * std::pow(size, _radius_exponent); };

  for (unsigned int n = 1; n <= _layout.switchSize(); ++n)
  {
    const unsigned int row = n - 1;
    _n[row] = n;
    _width[row] = 1.0;
    _concentration[row] = state(_layout.discreteComponent(n));
    _radius[row] = radius(static_cast<Real>(n));
  }
  unsigned int row = _layout.switchSize();
  for (const auto & cell : _layout.cells())
  {
    _n[row] = cell.center;
    _width[row] = cell.width;
    _concentration[row] = state(cell.component);
    _radius[row] = radius(cell.center);
    ++row;
  }
}

void
ClusterContinuumDistribution::threadJoin(const UserObject & /*y*/)
{
}
