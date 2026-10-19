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

#include "ClusterGroupingLayout.h"
#include "NodalVectorPostprocessor.h"

/**
 * Reports the full per-size cluster concentration distribution (one row per physical cluster size
 * n, from 1 to num_cluster_sizes) from a grouped L0/L1 cluster dynamics array variable, by
 * evaluating the same ClusterGroupingLayout::concentration() reconstruction the grouped kernels
 * use internally. Explicit sizes are read directly; grouped (large-n) sizes are reconstructed from
 * their bin's L0/L1 linear fit, so adjacent reconstructed rows within one bin lie on that bin's
 * single line rather than being independent data.
 */
class ClusterSizeDistributionGrouped : public NodalVectorPostprocessor
{
public:
  static InputParameters validParams();
  ClusterSizeDistributionGrouped(const InputParameters & parameters);

  void initialize() override;
  void execute() override;
  void threadJoin(const UserObject & y) override;

protected:
  const ArrayVariableValue & _clusters;
  const Real _r1;
  const Real _radius_exponent;
  const ClusterGroupingLayout _layout;
  const bool _enforce_group_nonnegative;
  const Real _group_nonnegative_tolerance_factor;

  VectorPostprocessorValue & _n;
  VectorPostprocessorValue & _concentration;
  VectorPostprocessorValue & _radius;
};
