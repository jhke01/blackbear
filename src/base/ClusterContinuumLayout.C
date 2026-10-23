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

#include "ClusterContinuumLayout.h"

#include "MooseError.h"

#include <algorithm>

ClusterContinuumLayout::ClusterContinuumLayout(unsigned int num_cluster_sizes,
                                               unsigned int switch_size,
                                               Real initial_width,
                                               Real width_ratio)
{
  if (num_cluster_sizes < 2)
    mooseError("ClusterContinuumLayout requires num_cluster_sizes >= 2.");
  if (switch_size < 2)
    mooseError("ClusterContinuumLayout requires switch_size >= 2.");
  if (switch_size >= num_cluster_sizes)
    mooseError("ClusterContinuumLayout requires switch_size < num_cluster_sizes; use a discrete "
               "layout when every size is an unknown.");
  if (initial_width <= 0.0)
    mooseError("ClusterContinuumLayout requires initial_width > 0.");
  if (width_ratio < 1.0)
    mooseError("ClusterContinuumLayout requires width_ratio >= 1.");

  _num_cluster_sizes = num_cluster_sizes;
  _switch_size = switch_size;

  // Component 0 is the monomer and components 1 .. switch_size - 1 are sizes 2 .. switch_size
  _component_count = switch_size;

  const Real domain_right = static_cast<Real>(num_cluster_sizes) + 0.5;
  Real left = static_cast<Real>(switch_size) + 0.5;
  Real width = initial_width;
  while (left < domain_right)
  {
    const Real remaining = domain_right - left;
    Real w = std::min(width, remaining);
    // Merge a leftover sliver of less than half a cell into this cell
    if (remaining > w && remaining - w < 0.5 * w)
      w = remaining;

    ClusterContinuumCell cell;
    cell.left = left;
    cell.right = left + w;
    cell.center = 0.5 * (cell.left + cell.right);
    cell.width = w;
    cell.first_moment = 0.5 * (cell.right * cell.right - cell.left * cell.left);
    cell.component = _component_count++;
    _cells.push_back(cell);

    left = cell.right;
    width *= width_ratio;
  }
}

unsigned int
ClusterContinuumLayout::discreteComponent(unsigned int n) const
{
  if (!isDiscreteSize(n))
    mooseError("Requested discrete component for cluster size ", n, ", which is not discrete.");
  return n - 1;
}

unsigned int
ClusterContinuumLayout::cellForSize(Real x) const
{
  if (_cells.empty())
    mooseError("ClusterContinuumLayout has no cells.");
  if (x <= _cells.front().left)
    return 0;
  if (x >= _cells.back().right)
    return numCells() - 1;

  // Cell right edges increase, so the first cell whose right edge exceeds x contains it
  const auto it = std::upper_bound(_cells.begin(),
                                   _cells.end(),
                                   x,
                                   [](Real value, const ClusterContinuumCell & c)
                                   { return value < c.right; });
  return static_cast<unsigned int>(it - _cells.begin());
}
