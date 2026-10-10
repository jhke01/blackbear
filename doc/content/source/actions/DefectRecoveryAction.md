# DefectRecoveryAction

!syntax description /DefectDynamics/Recovery/DefectRecoveryAction

## Description

The `DefectRecoveryAction` is a MOOSE action that describes the shared dislocation network
parameters (initial density, recovery rate constant, Burgers vector) used by every species'
[GenericDislocationRecoveryNodalKernel.md] in the enclosing `DefectDynamics` block, when this
optional sub-block is given. This is invoked using the `Recovery` sub-block within
[DefectDynamics](syntax/DefectDynamics/index.md). See that page for a description, example
usage, and parameters.
