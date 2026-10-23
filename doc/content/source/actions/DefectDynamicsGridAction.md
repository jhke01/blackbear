# DefectDynamicsGridAction

!syntax description /DefectDynamics/Grid/DefectDynamicsGridAction

## Description

The `DefectDynamicsGridAction` is a MOOSE action that describes a grouped-tail layout shared by
every species in the enclosing `DefectDynamics` block. It is deprecated: each species now sets
its own grouping, with `use_grouping`, `explicit_max`, `group_spacing`, `group_linear_width`
and `group_geometric_ratio` on its own sub-block. A `Grid` sub-block still works and supplies the
grouping for any species that did not set its own, and it prints a deprecation warning. See
[DefectDynamics](syntax/DefectDynamics/index.md) for the current inputs.
