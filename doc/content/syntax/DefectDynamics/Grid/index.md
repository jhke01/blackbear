# DefectDynamics/Grid

The `[Grid]` sub-block of [DefectDynamics](syntax/DefectDynamics/index.md) is deprecated. It describes a grouped-tail layout that every species of the enclosing block
shares. Each species now sets its own grouping (`use_grouping`, `explicit_max`, `group_spacing`, `group_linear_width` and `group_geometric_ratio`) on its own sub-block. A
`[Grid]` sub-block still works: it supplies the grouping for any species that did not set its own, and it prints a deprecation warning. See [DefectDynamicsGridAction.md]
for the parameters.

## Constructed Objects

This block holds only parameters. [DefectDynamics](syntax/DefectDynamics/index.md) reads them when it builds the species.

## Example Input Syntax

!syntax parameters /DefectDynamics/Grid/DefectDynamicsGridAction
