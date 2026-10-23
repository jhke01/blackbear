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

#include "ClusterContinuumMoment.h"

#include <cmath>

registerMooseObject("BlackBearApp", ClusterContinuumMoment);

namespace
{
/// Integral of x^exponent from left to right
Real
powerIntegral(Real left, Real right, Real exponent)
{
  const Real p1 = exponent + 1.0;
  return (std::pow(right, p1) - std::pow(left, p1)) / p1;
}
} // namespace

InputParameters
ClusterContinuumMoment::validParams()
{
  InputParameters params = NodalPostprocessor::validParams();
  params.set<bool>("unique_node_execute") = true;
  params.addClassDescription(
      "Computes the number density, the point-defect content or the mean radius of a cluster size "
      "distribution that is discrete up to switch_size and a Fokker-Planck continuum above it.");
  params.addRequiredCoupledVar("clusters", "Cluster array concentration variable.");
  MooseEnum quantity("density content average_radius");
  params.addRequiredParam<MooseEnum>(
      "quantity",
      quantity,
      "density: sum of C_n over V_at [1/m^3]; content: sum of n C_n over V_at [1/m^3]; "
      "average_radius: number-weighted mean radius [m] (r_n = r1 n^radius_exponent).");
  params.addRequiredParam<unsigned int>("num_cluster_sizes",
                                        "Largest cluster size represented, N.");
  params.addRequiredParam<unsigned int>(
      "switch_size", "Largest cluster size that is an individual unknown; larger sizes are cells.");
  params.addParam<Real>("continuum_initial_width", 1.0, "Width of the first continuum cell.");
  params.addParam<Real>(
      "continuum_width_ratio", 1.04, "Growth factor of the cell width from one cell to the next.");
  params.addRangeCheckedParam<unsigned int>(
      "n_minimum",
      1,
      "n_minimum > 0",
      "Smallest cluster size included (1 includes the monomer).");
  params.addParam<Real>("atomic_volume", 1.0, "Atomic volume [m^3] for density and content.");
  params.addParam<Real>("r1", 1.0, "Monomer radius [m] in r_n = r1*n^radius_exponent.");
  params.addParam<Real>("radius_exponent",
                        0.5,
                        "Exponent in r_n = r1*n^radius_exponent: 1/2 for loops, 1/3 for voids.");
  return params;
}

ClusterContinuumMoment::ClusterContinuumMoment(const InputParameters & parameters)
  : NodalPostprocessor(parameters),
    _clusters(coupledArrayValue("clusters")),
    _quantity(getParam<MooseEnum>("quantity").getEnum<Quantity>()),
    _r1(getParam<Real>("r1")),
    _radius_exponent(getParam<Real>("radius_exponent")),
    _n_minimum(getParam<unsigned int>("n_minimum")),
    _atomic_volume(getParam<Real>("atomic_volume")),
    _layout(getParam<unsigned int>("num_cluster_sizes"),
            getParam<unsigned int>("switch_size"),
            getParam<Real>("continuum_initial_width"),
            getParam<Real>("continuum_width_ratio")),
    _number(0.0),
    _content(0.0),
    _radius_weighted(0.0)
{
  if (_atomic_volume <= 0.0)
    paramError("atomic_volume", "must be positive.");
  if (_r1 <= 0.0 || _radius_exponent <= 0.0)
    mooseError(type(), " requires r1 > 0 and radius_exponent > 0.");
  if (getArrayVar("clusters", 0)->count() != _layout.componentCount())
    paramError("clusters",
               "The continuum layout requires ",
               _layout.componentCount(),
               " array components, but the coupled variable has ",
               getArrayVar("clusters", 0)->count(),
               ". The layout parameters must match those of the cluster variable.");
}

void
ClusterContinuumMoment::initialize()
{
  _number = 0.0;
  _content = 0.0;
  _radius_weighted = 0.0;
}

void
ClusterContinuumMoment::execute()
{
  const RealEigenVector & state = _clusters[_qp];
  const unsigned int last_discrete = std::min(_layout.switchSize(), _layout.numClusterSizes());
  for (unsigned int n = _n_minimum; n <= last_discrete; ++n)
  {
    // Concentrations enter as max(C, 0), as in the kernel
    const Real c = std::max(state(_layout.discreteComponent(n)), 0.0);
    _number += c;
    _content += static_cast<Real>(n) * c;
    _radius_weighted += _r1 * std::pow(static_cast<Real>(n), _radius_exponent) * c;
  }

  // The lower edge of size n_minimum; the cells begin at switch_size + 1/2
  const Real threshold_edge = static_cast<Real>(_n_minimum) - 0.5;
  for (const auto & cell : _layout.cells())
  {
    const Real left = std::max(cell.left, threshold_edge);
    if (cell.right <= left)
      continue;
    const Real c = std::max(state(cell.component), 0.0);
    _number += c * (cell.right - left);
    _content += c * powerIntegral(left, cell.right, 1.0);
    _radius_weighted += c * _r1 * powerIntegral(left, cell.right, _radius_exponent);
  }
}

void
ClusterContinuumMoment::finalize()
{
  gatherSum(_number);
  gatherSum(_content);
  gatherSum(_radius_weighted);
}

void
ClusterContinuumMoment::threadJoin(const UserObject & y)
{
  const auto & rhs = static_cast<const ClusterContinuumMoment &>(y);
  _number += rhs._number;
  _content += rhs._content;
  _radius_weighted += rhs._radius_weighted;
}

Real
ClusterContinuumMoment::getValue() const
{
  switch (_quantity)
  {
    case Quantity::DENSITY:
      return _number / _atomic_volume;
    case Quantity::CONTENT:
      return _content / _atomic_volume;
    case Quantity::AVERAGE_RADIUS:
      return _number > 0.0 ? _radius_weighted / _number : 0.0;
  }
  return 0.0;
}
