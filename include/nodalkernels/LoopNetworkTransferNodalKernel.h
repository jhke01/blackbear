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

#include "ArrayNodalKernel.h"

/**
 * Loss of the interstitial loops of a cluster size distribution to the dislocation network, the
 * counterpart of the line-length source of CellularNetworkDislocationDensity when that postprocessor
 * is given a line_length_postprocessor.
 *
 * A loop climbs with the velocity v_l = 2 Z_iv D_i c_i / b (Rest 2006, Eq 16), where c_i is the
 * monomer (array component 0) of this same variable, and reaches the network after the distance
 * L = C_A C_rho sqrt(pi / (rho f(nu))) between its walls (Rest Eq 10). Every loop of at least
 * minimum_size therefore leaves its size at the rate
 *
 *   k = v_l / L = v_l sqrt(f(nu)/pi) / (C_A C_rho) sqrt(rho),
 *
 * independent of its size, so the residual is k C_n on every component from minimum_size up (discrete
 * sizes and continuum cells alike). The line length 2 pi r_n C_n / V_at of the loops that leave is
 * added to the network density by the postprocessor, so the total line length is conserved by the
 * pair. The density rho is read from a postprocessor and is constant within a step; the dependence of
 * k on the monomer c_i is part of the Jacobian.
 */
class LoopNetworkTransferNodalKernel : public ArrayNodalKernel
{
public:
  static InputParameters validParams();
  LoopNetworkTransferNodalKernel(const InputParameters & parameters);

protected:
  virtual void computeQpResidual(RealEigenVector & residual) override;
  virtual void computeQpJacobian() override;

  /// Rate k [1/s] at which one loop joins the network, for the monomer concentration c
  Real rate(Real c) const;

  const Real _coefficient; ///< 2 Z_iv D_i sqrt(f/pi) / (b C_A C_rho)
  const unsigned int _minimum_size;
  const PostprocessorValue & _rho;
};
