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
#include "NodalPostprocessor.h"

/**
 * Moments of a cluster size distribution that is discrete up to switch_size and a Fokker-Planck
 * continuum above it (ClusterContinuumLayout), restricted to sizes of at least n_minimum:
 * the number density sum C_n / V_at, the point-defect content sum n C_n / V_at, or the
 * number-weighted mean radius sum r_n C_n / sum C_n with r_n = r1 n^radius_exponent.
 *
 * A continuum cell holds a concentration per unit size, so its moment is the integral over the
 * part of the cell above n_minimum - 1/2 (the lower edge of size n_minimum), as in the SUNDIALS
 * Pokor benchmark. The monomer is size 1.
 */
class ClusterContinuumMoment : public NodalPostprocessor
{
public:
  static InputParameters validParams();
  ClusterContinuumMoment(const InputParameters & parameters);

  void initialize() override;
  void execute() override;
  void finalize() override;
  void threadJoin(const UserObject & y) override;
  Real getValue() const override;

protected:
  enum class Quantity
  {
    DENSITY,
    CONTENT,
    AVERAGE_RADIUS
  };

  const ArrayVariableValue & _clusters;
  const Quantity _quantity;
  const Real _r1;
  const Real _radius_exponent;
  const unsigned int _n_minimum;
  const Real _atomic_volume;
  const ClusterContinuumLayout _layout;

  /// sum of C_n (the number of clusters), sum of n C_n, and sum of r_n C_n
  Real _number;
  Real _content;
  Real _radius_weighted;
};
