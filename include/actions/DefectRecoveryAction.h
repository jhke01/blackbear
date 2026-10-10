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
 * Holds no logic of its own; it exists only so DefectDynamicsAction can look it up
 * (ActionWarehouse::getActions<DefectRecoveryAction>()) and read its validated parameters while
 * building the DislocationRecoveryNodalKernel instances.
 */
class DefectRecoveryAction : public Action
{
public:
  static InputParameters validParams();

  DefectRecoveryAction(const InputParameters & params);

  virtual void act() override {}
};
