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

#include "DefectRecoveryAction.h"

registerMooseAction("BlackBearApp", DefectRecoveryAction, "meta_action");

InputParameters
DefectRecoveryAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription(
      "The shared dislocation network parameters (initial density, recovery rate constant, "
      "Burgers vector) used by every species' DislocationRecoveryNodalKernel in the enclosing "
      "DefectDynamics block, when this optional sub-block is given.");
  params.addRequiredParam<Real>(
      "rho0", "Initial dislocation density rho0 [m^-2] of the shared dislocation network.");
  params.addRequiredParam<Real>(
      "K",
      "Recovery rate constant K [m^-1*s^-1] in d(rho)/dt = -K*b^2*rho^(3/2) (Pokor et al. 2004, "
      "Eq 8). Note: the paper's own text states K in bare m^-1, which is dimensionally "
      "inconsistent with this governing equation; treat K as a fit parameter in self-consistent "
      "SI units, not a literal copy of the paper's quoted numbers.");
  params.addRequiredParam<Real>(
      "burgers_vector",
      "Burgers vector b [m] of the dislocation network itself -- distinct from any loop "
      "species' own geometric burgers_vector parameter.");
  return params;
}

DefectRecoveryAction::DefectRecoveryAction(const InputParameters & params) : Action(params)
{
}
