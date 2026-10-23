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

#include "Moose.h"

#include <vector>

/**
 * One finite-volume cell of the continuum tail. The cell covers the real-valued size interval
 * [left, right] and its unknown is the cell-average concentration per unit size.
 */
struct ClusterContinuumCell
{
  Real left = 0.0;
  Real right = 0.0;
  Real center = 0.0;
  Real width = 0.0;

  /// Integral of n dn over the cell, 0.5*(right^2 - left^2): the number of point defects per unit
  /// concentration held by the cell
  Real first_moment = 0.0;

  /// Array component that holds the cell average
  unsigned int component = libMesh::invalid_uint;
};

/**
 * Array layout for a cluster size distribution that is discrete up to switch_size and continuous
 * above it, the alternative to ClusterGroupingLayout's L0/L1 bins. The continuum equation is
 * discretized by finite volumes (Pokor et al. 2004, Section 4.1).
 *
 * Physical cluster size n is one-based, as in ClusterGroupingLayout. Array component layout:
 *   component 0 -> monomer
 *   components 1 .. switch_size-1 -> sizes 2 .. switch_size, one component each
 *   components switch_size .. switch_size+cells-1 -> one component per continuum cell
 *
 * The cells start at switch_size + 0.5 with width initial_width, the width grows by width_ratio
 * from one cell to the next, and the last cell ends at num_cluster_sizes + 0.5. When the last
 * cell would be left shorter than half of the previous width, it is merged into the previous
 * cell. These are the rules of the SUNDIALS Pokor benchmark (Layout::build in
 * pokor_cluster_dynamics.cpp), so the two codes have the same number of unknowns.
 */
class ClusterContinuumLayout
{
public:
  ClusterContinuumLayout() = default;

  ClusterContinuumLayout(unsigned int num_cluster_sizes,
                         unsigned int switch_size,
                         Real initial_width,
                         Real width_ratio);

  unsigned int numClusterSizes() const { return _num_cluster_sizes; }
  unsigned int switchSize() const { return _switch_size; }
  unsigned int componentCount() const { return _component_count; }
  unsigned int monomerComponent() const { return 0; }
  unsigned int numCells() const { return static_cast<unsigned int>(_cells.size()); }
  const std::vector<ClusterContinuumCell> & cells() const { return _cells; }
  const ClusterContinuumCell & cell(unsigned int index) const { return _cells[index]; }

  bool containsSize(unsigned int n) const { return n >= 1 && n <= _num_cluster_sizes; }
  bool isDiscreteSize(unsigned int n) const { return containsSize(n) && n <= _switch_size; }

  /// Array component of a discrete size 1 <= n <= switch_size
  unsigned int discreteComponent(unsigned int n) const;

  /// Index of the cell that contains the (real-valued) size x, between 0 and numCells() - 1
  unsigned int cellForSize(Real x) const;

private:
  unsigned int _num_cluster_sizes = 0;
  unsigned int _switch_size = 0;
  unsigned int _component_count = 1;
  std::vector<ClusterContinuumCell> _cells;
};
