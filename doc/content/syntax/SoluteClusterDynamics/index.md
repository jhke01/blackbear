# SoluteClusterDynamics

!syntax description /SoluteClusterDynamics/SoluteClusterDynamicsAction

## Description

The `SoluteClusterDynamics` block is a front end to the solute-clustering kernel,
[GenericClusterDynamicsGroupedNodalKernel.md] (commonly used for Cu precipitation in
reactor-pressure-vessel steels). It builds the solute species' grouped array variable, initial
condition, time-derivative kernel, dynamics kernel, and monomer-concentration/total-density
postprocessors from the physical parameters given below, instead of requiring
`[GroupedVariables]`/`[NodalKernels]` to be hand-typed. Unlike [DefectDynamics/index.md], the Cu
kernel already derives its own rate coefficients internally from the raw physical parameters, so
this action is mostly a thin pass-through of the kernel's own parameters, with one piece of
genuine added value: automatic point-defect coupling.

## Automatic Point-Defect Coupling

When a `[DefectDynamics]` block is ALSO present in the same input, `SoluteClusterDynamics`
automatically looks for that block's vacancy species -- identified by `geometry = void`, falling
back to a species literally named `vacancy` if no `geometry = void` species is present (covering a
Pokor-style setup where vacancy is modeled as `geometry = loop`, see [DefectDynamics/index.md]) --
and, when found, wires
[!param](/NodalKernels/ClusterDynamicsGroupedNodalKernel/vacancy_coupled_var) and
[!param](/NodalKernels/ClusterDynamicsGroupedNodalKernel/vacancy_formation_energy_eV) onto the Cu
kernel it builds. No separate switch exists to turn this on: declaring both blocks together in
the same input **is** how a user studying solute clusters activates concurrent point-defect
clustering and its effect on solute diffusion. A `[SoluteClusterDynamics]` block with no
`[DefectDynamics]` block present builds the solute species exactly as it would otherwise, with
the coupling simply inactive.

```
[DefectDynamics]
  temperature = 600
  atomic_volume = 1.1782924e-29
  lattice_parameter = 0.36e-9
  dpa_rate = 1.5e-4
  cascade_efficiency = 0.1
  dislocation_density = 1e15

  [vacancy]
    geometry = void
    D0 = 1.19073105525e-8
    Em = 0.7644
    Ef = 1.6
    Eb2 = 0.5
    exponent = 0.67
    num_cluster_sizes = 30
    explicit_max = 6
    group_spacing = linear
    group_linear_width = 4
    group_geometric_ratio = 1.5
  []
[]

[SoluteClusterDynamics]
  initial_concentration = 0.003
  rate_model = interfacial_energy
  diffusivity_model = constant
  monomer_diffusivity = 1.0e-21
  temperature = 600
  atomic_volume = 1.1782924e-29
  sigma = 0.37
  Omega_kB_K = 6255
  DeltaS_kB = 0.866

  num_cluster_sizes = 50000
  explicit_max = 200
  group_spacing = geometric
  group_geometric_ratio = 1.1
[]
```

Here, `vacancy_coupled_var` is set to the vacancy species' own array variable
(`vacancy_clusters`, following [DefectDynamics/index.md]'s own naming convention) and
`vacancy_formation_energy_eV` to `1.6` (the vacancy species' own `Ef`), with neither typed
directly on the `[SoluteClusterDynamics]` block.

If more than one `[DefectDynamics]` block is present, or if the vacancy species cannot be
unambiguously identified (for example, two `geometry = void` species in the same
`[DefectDynamics]` block), `SoluteClusterDynamics` raises a clear error rather than silently
picking one.

## Grouping and Resolvable Clusters

The solute array has its own grouping, set on the `[SoluteClusterDynamics]` block with the same
names a `[DefectDynamics]` species uses:
[!param](/SoluteClusterDynamics/SoluteClusterDynamicsAction/use_grouping),
[!param](/SoluteClusterDynamics/SoluteClusterDynamicsAction/explicit_max),
[!param](/SoluteClusterDynamics/SoluteClusterDynamicsAction/group_spacing),
[!param](/SoluteClusterDynamics/SoluteClusterDynamicsAction/group_linear_width) and
[!param](/SoluteClusterDynamics/SoluteClusterDynamicsAction/group_geometric_ratio). It does not
depend on the grouping of any `[DefectDynamics]` species.

[!param](/SoluteClusterDynamics/SoluteClusterDynamicsAction/resolvable_radius) works as it does on
a `[DefectDynamics]` species. It gives the smallest cluster radius [m] that a characterization
technique can see, and adds the outputs `<name>_resolvable_density` and
`<name>_resolvable_avg_radius`, which count only clusters at or above that radius. It needs
[!param](/SoluteClusterDynamics/SoluteClusterDynamicsAction/atomic_volume) to convert the radius to
a cluster size. It replaces the older `density_n_minimum`, which still works but is deprecated.

## Fokker-Planck Continuum Tail

Instead of grouping the large sizes in bins, the solute array can store them as a Fokker-Planck
continuum (Pokor et al. 2004, Section 4.1) by setting
[!param](/SoluteClusterDynamics/SoluteClusterDynamicsAction/tail) `= fokker_planck`. Sizes up to
[!param](/SoluteClusterDynamics/SoluteClusterDynamicsAction/switch_size) are individual unknowns, and larger
sizes are finite-volume cells that start at `switch_size` + 1/2 with the width
[!param](/SoluteClusterDynamics/SoluteClusterDynamicsAction/continuum_initial_width) and grow by
[!param](/SoluteClusterDynamics/SoluteClusterDynamicsAction/continuum_width_ratio) from one cell to the next.
The block then builds a [ContinuumClusterDynamicsNodalKernel.md] in place of the grouped kernel, the monomer
postprocessor reads the array's first component, and the densities and the mean radius come from
[ClusterContinuumMoment.md]. The automatic coupling to the vacancies of a `[DefectDynamics]` block works as
for the grouped array.

The grouping parameters (`use_grouping`, `explicit_max`, `group_*`) are not allowed with the continuum, and the
binding-energy tables are not supported (`interfacial_energy` and `capillary` are). The block warns when
`switch_size` is below the size at which clusters become resolvable
([!param](/SoluteClusterDynamics/SoluteClusterDynamicsAction/resolvable_radius)), because the cut of the
resolvable population then falls inside a continuum cell.

The continuum is useful where the memory of the grouped array matters: the cost of the Jacobian storage grows
with the square of the number of unknowns of all the arrays of the problem, and the continuum needs far fewer
unknowns for the same accuracy of the point-defect populations of [DefectDynamics/index.md]. For the Cu
precipitation benchmark with 100,000 sizes, which the grouped array already represents with 338 unknowns, the
continuum (311 to 584 unknowns) needs about the same memory (0.03 to 0.04 GB). The two schemes are equally
accurate once the time step is controlled: with the step limited to $0.005\,t + 10^{-4}$ s, both agree with the
SUNDIALS reference to within 0.4 percent (grouped) and 1.3 percent (continuum) in the density, and 0.2 and 0.4
percent in the radius, for 3095 steps; the continuum is about 6.5 times faster (21 s against 140 s). With the
default `IterationAdaptiveDT` alone, both are off by a factor of 3 in the density at 0.01 s and by 18 percent at
1 s, because the steps are too large for the early rise of the cluster density.

## Time Step for the Onset of the Resolvable Population

With a `resolvable_radius`, the resolvable density rises by many orders of magnitude in a short interval, and a time step that is only capped in proportion to
the time cannot follow it (see [DefectDynamics/index.md]). Setting
[!param](/SoluteClusterDynamics/SoluteClusterDynamicsAction/timestep_relative_change) adds the postprocessor
`<block name>_timestep_limit` ([RelativeChangeTimestepLimit.md]), which limits the relative change of the resolvable density to that fraction per step,
bounded above by [!param](/SoluteClusterDynamics/SoluteClusterDynamicsAction/timestep_cap_fraction) $\times\, t\,+$
[!param](/SoluteClusterDynamics/SoluteClusterDynamicsAction/timestep_cap_minimum). The time stepper then names it as its
`timestep_limiting_postprocessor`. It needs `resolvable_radius`.

## Example Input Syntax

In this example, no `[DefectDynamics]` block is present, so a solute-only study is built with the
coupling simply inactive.

!listing test/tests/solute_cluster_dynamics/solute_only_action.i block=SoluteClusterDynamics

!syntax parameters /SoluteClusterDynamics/SoluteClusterDynamicsAction

## See Also

- [GenericClusterDynamicsGroupedNodalKernel.md]
- [ContinuumClusterDynamicsNodalKernel.md]
- [ClusterContinuumMoment.md]
- [DefectDynamics/index.md]
