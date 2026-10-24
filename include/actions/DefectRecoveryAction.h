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

#include "Action.h"

/**
 * The [DefectDynamics/Recovery] sub-block: the shared dislocation network's initial density
 * (rho0), recovery rate constant (K), and Burgers vector (burgers_vector), used identically by
 * every species present when this block is given (Pokor et al. 2004, Eqs 7-8 use one
 * un-subscripted dislocation density for both species' sink terms). Optional -- zero or one
 * instance; when absent, every species keeps its existing constant dislocation sink.
 *
 * With model = cellular_network the density is instead the cellular network dislocation density of
 * Rest (2006), advanced by a CellularNetworkDislocationDensity postprocessor that
 * DefectDynamicsAction builds from the interstitial loop species' monomer concentration; K is then
 * not used.
 *
 * Holds no logic of its own beyond validating its parameters; it exists so DefectDynamicsAction can
 * look it up (ActionWarehouse::getActions<DefectRecoveryAction>()) and read them while building the
 * DislocationRecoveryNodalKernel instances.
 */
class DefectRecoveryAction : public Action
{
public:
  static InputParameters validParams();

  DefectRecoveryAction(const InputParameters & params);

  virtual void act() override {}

  /// True when the density is the cellular network model of Rest (2006)
  bool cellularNetwork() const;
};
