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
  params.addParam<Real>(
      "K",
      "Recovery rate constant K [m^-1*s^-1] in d(rho)/dt = -K*b^2*rho^(3/2) (Pokor et al. 2004, "
      "Eq 8). Required for model = pokor, not allowed for model = cellular_network. Note: the paper's own text states K in bare m^-1, which is dimensionally "
      "inconsistent with this governing equation; treat K as a fit parameter in self-consistent "
      "SI units, not a literal copy of the paper's quoted numbers.");
  params.addRequiredParam<Real>(
      "burgers_vector",
      "Burgers vector b [m] of the dislocation network itself -- distinct from any loop "
      "species' own geometric burgers_vector parameter.");
  MooseEnum model("pokor cellular_network", "pokor");
  params.addParam<MooseEnum>(
      "model",
      model,
      "pokor: the density recovers by climb, d(rho)/dt = -K b^2 rho^(3/2) (Pokor et al. 2004), "
      "in closed form. cellular_network: the density is created by interstitial loops that join "
      "the network and lost by self-annihilation (Rest 2006, Eqs 14-18), advanced from the "
      "interstitial monomer concentration of the loop species.");
  params.addParam<Real>("bias",
                        "Relative interstitial/vacancy bias Z_iv of Rest (2006) Eq 16, "
                        "model = cellular_network only. Material dependent (5e-4 for UO2); not "
                        "known for steel.");
  params.addParam<Real>("ca_crho",
                        2.0,
                        "Product C_A C_rho of the cell geometry (Rest 2006 Eqs 10, 14), "
                        "model = cellular_network only. C_A = 3 for cubic cells, C_rho of order one.");
  params.addParam<Real>("poisson_ratio",
                        0.3,
                        "Poisson ratio nu of f(nu) = (1 - nu/2)/(1 - nu), model = cellular_network "
                        "only.");
  params.addParam<Real>("grain_size",
                        0.0,
                        "Grain size d_g [m] of the grain-boundary loss 4 v_l rho/d_g, "
                        "model = cellular_network only. 0 turns the term off.");
  params.addParam<bool>(
      "loop_transfer",
      false,
      "model = cellular_network only. If true, the loops of the interstitial species join the network "
      "at the rate k = v_l / L (L is the cell size of Rest's Eq 10) and carry their line length into "
      "the network density: LoopNetworkTransferNodalKernel removes them from the size distribution "
      "and the network density gains their line length. This removes the double counting of line "
      "length between the explicit loops and the network. If false, the network is fed by the "
      "steady-state loop density of Rest's Eq 18 and the explicit loops are not affected.");
  params.addRangeCheckedParam<unsigned int>(
      "transfer_minimum_size",
      2,
      "transfer_minimum_size >= 2",
      "loop_transfer only: smallest loop size that can join the network.");
  params.addParam<std::string>(
      "interstitial_species",
      "Name of the species sub-block whose monomer concentration drives the network, "
      "model = cellular_network only. Needed when more than one species has geometry = loop.");
  return params;
}

DefectRecoveryAction::DefectRecoveryAction(const InputParameters & params) : Action(params)
{
  if (cellularNetwork())
  {
    if (isParamValid("K"))
      paramError("K", "is not used with model = cellular_network.");
    if (!isParamValid("bias"))
      paramError("bias", "is required with model = cellular_network.");
    if (getParam<Real>("bias") <= 0.0)
      paramError("bias", "must be positive.");
    if (getParam<Real>("ca_crho") <= 0.0)
      paramError("ca_crho", "must be positive.");
    if (getParam<Real>("poisson_ratio") < 0.0 || getParam<Real>("poisson_ratio") >= 1.0)
      paramError("poisson_ratio", "must be in [0, 1).");
    if (getParam<Real>("grain_size") < 0.0)
      paramError("grain_size", "must not be negative.");
    if (!getParam<bool>("loop_transfer") && parameters().isParamSetByUser("transfer_minimum_size"))
      paramError("transfer_minimum_size", "is only used with loop_transfer = true.");
  }
  else
  {
    if (!isParamValid("K"))
      paramError("K", "is required with model = pokor.");
    for (const auto & name :
         {"bias", "interstitial_species", "loop_transfer", "transfer_minimum_size"})
      if (parameters().isParamSetByUser(name))
        paramError(name, "is only used with model = cellular_network.");
  }
}

bool
DefectRecoveryAction::cellularNetwork() const
{
  return getParam<MooseEnum>("model") == "cellular_network";
}
