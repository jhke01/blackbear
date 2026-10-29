# DefectDynamics/Recovery

The optional `[Recovery]` sub-block of [DefectDynamics](syntax/DefectDynamics/index.md) describes the dislocation network that is shared by the
species of the enclosing block. It replaces the constant dislocation sink of each species with a dislocation density that changes in time. The
[!param](/DefectDynamics/Recovery/DefectRecoveryAction/model) parameter selects Pokor's recovery law (`pokor`, the default) or the cellular
network model of [!cite](rest2006) (`cellular_network`). See [DefectRecoveryAction.md] for the parameters and the models, and
[theory/defect_dynamics.md] for the equations.

## Constructed Objects

This block holds only parameters. [DefectDynamics](syntax/DefectDynamics/index.md) reads them and adds one [GenericDislocationRecoveryNodalKernel.md]
for each species and, with `model = cellular_network`, a [CellularNetworkDislocationDensity.md] postprocessor, and with `loop_transfer = true` a
[LoopNetworkTransferNodalKernel.md].

## Example Input Syntax

!listing test/tests/defect_dynamics/defect_dynamics_cellular_network_action.i block=DefectDynamics/Recovery

!syntax parameters /DefectDynamics/Recovery/DefectRecoveryAction
