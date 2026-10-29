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

#include "GeneralPostprocessor.h"

/**
 * Cellular network dislocation density rho_N(t) [m^-2] of Rest (J. Nucl. Mater. 349 (2006) 150),
 * advanced at the end of each time step. The interstitial loop density follows from its steady
 * state (Rest Eq 18, n_l = k sqrt(rho_N) with k^2 = sqrt(2) b c_i / (4 pi Z_iv Omega^(5/3))), so
 * Rest Eq 14 becomes, with u = sqrt(rho_N), v_l = 2 Z_iv D_i c_i / b (Eq 16),
 *
 *   du/dt = (pi v_l k)/2 - (2 v_l/d_g) u - (C v_l/2) u^2,  C = (f(nu)/pi)^(1/2)/(C_A C_rho),
 *
 * a Riccati equation that is solved exactly over the step for a constant c_i (the mean of its
 * values at the two ends of the step). The step is therefore unconditionally stable and keeps
 * rho_N positive; the accuracy is that of holding c_i constant, which is controlled by the time step
 * (RelativeChangeTimestepLimit on this postprocessor).
 *
 * The value is the density at the end of the last completed step. Kernels that read it
 * (DislocationRecoveryNodalKernel with dislocation_density_postprocessor) see it as a constant
 * within a step.
 */
class CellularNetworkDislocationDensity : public GeneralPostprocessor
{
public:
  static InputParameters validParams();
  CellularNetworkDislocationDensity(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override;
  virtual void finalize() override {}
  virtual Real getValue() const override { return _rho; }

  /// rho after a time dt with the interstitial site fraction c held constant, from rho0
  static Real advance(Real rho0,
                      Real dt,
                      Real c_i,
                      Real diffusivity,
                      Real burgers_vector,
                      Real atomic_volume,
                      Real z_iv,
                      Real ca_crho,
                      Real poisson_ratio,
                      Real grain_size);

protected:
  const PostprocessorValue & _c_i;
  const Real _D0;
  const Real _Em_eV;
  const Real _temperature;
  const Real _burgers_vector;
  const Real _atomic_volume;
  const Real _rho0;
  const Real _z_iv;
  const Real _ca_crho;
  const Real _poisson_ratio;
  const Real _grain_size;
  const Real _diffusivity;

  /// State: the density and the interstitial site fraction at the end of the last step
  Real & _rho;
  Real & _c_i_old;
};
