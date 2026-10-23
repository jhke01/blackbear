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

#include "ClusterContinuumLayout.h"
#include "NodalVectorPostprocessor.h"

/**
 * Reports the size distribution of a ClusterContinuumLayout variable, one row per discrete size
 * and one row per continuum cell: the size (the cell center for a cell), the width of the size
 * interval the row covers (1 for a discrete size), the concentration per unit size, and the
 * radius at that size. Multiplying a row's concentration by its width gives the number of
 * clusters in the interval.
 */
class ClusterContinuumDistribution : public NodalVectorPostprocessor
{
public:
  static InputParameters validParams();
  ClusterContinuumDistribution(const InputParameters & parameters);

  void initialize() override;
  void execute() override;
  void threadJoin(const UserObject & y) override;

protected:
  const ArrayVariableValue & _clusters;
  const Real _r1;
  const Real _radius_exponent;
  const ClusterContinuumLayout _layout;

  VectorPostprocessorValue & _n;
  VectorPostprocessorValue & _width;
  VectorPostprocessorValue & _concentration;
  VectorPostprocessorValue & _radius;
};
