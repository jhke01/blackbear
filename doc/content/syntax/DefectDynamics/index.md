# DefectDynamics

!syntax description /DefectDynamics/DefectDynamicsAction

## Description

The `DefectDynamics` block is a species-organized front end to the grouped defect-dynamics
kernels ([GenericVacancyVoidGroupedNodalKernel.md], [GenericInterstitialLoopGroupedNodalKernel.md],
[PointDefectRecombinationNodalKernel.md], [GenericDislocationRecoveryNodalKernel.md],
[GenericCascadeSpectrumNodalKernel.md]; see [theory/defect_dynamics.md] for the underlying
equations). One or two named species sub-blocks (commonly `[vacancy]` and `[interstitial]`,
though the name is only a label) describe each species' own diffusivity, capillary
formation-energy law, (for a loop) dislocation bias, its own maximum cluster size and its own
grouping of the size array (explicit size, bin spacing, growth ratio; see below). One optional
`[Recovery]` sub-block describes a dislocation network shared by both species whose density
evolves in time instead of staying fixed. A species also
turns on cascade seeding of cluster sizes 2-4 simply by setting a nonzero cascade fraction on its
own sub-block.

Void and loop populations routinely grow to very different physical sizes under the same
irradiation conditions, so the maximum cluster size is set per species via
[!param](/DefectDynamics/DefectSpeciesAction/num_cluster_sizes) (the cap given directly, as a
point-defect count) or
[!param](/DefectDynamics/DefectSpeciesAction/max_cluster_radius_m) (the cap given as a physical
radius in meters, converted internally using that species' own monomer radius law); exactly one
of the two must be given per species.

### Grouping

Each species describes the grouping of its own array under its own sub-block, with the same
parameter names the grouped kernels use: [!param](/DefectDynamics/DefectSpeciesAction/use_grouping),
[!param](/DefectDynamics/DefectSpeciesAction/explicit_max),
[!param](/DefectDynamics/DefectSpeciesAction/group_spacing),
[!param](/DefectDynamics/DefectSpeciesAction/group_linear_width) and
[!param](/DefectDynamics/DefectSpeciesAction/group_geometric_ratio). A coarser grouping for the
vacancy clusters than for the loops is therefore just a different value under `[vacancy]`. The
opposite-species parameters of each kernel are filled in from the partner's own grouping.
Setting [!param](/DefectDynamics/DefectSpeciesAction/use_grouping) to false gives every cluster
size its own unknown; because the cost grows roughly with the square of the unknown count, a
species without grouping is limited to 3000 sizes.

!alert warning title=The `[Grid]` sub-block is deprecated
The former shared `[Grid]` sub-block still works and supplies the grouping for any species that
does not set its own, but it prints a deprecation warning. Move its parameters
(`spacing` becomes `group_spacing`, `linear_width` becomes `group_linear_width`,
`geometric_ratio` becomes `group_geometric_ratio`) into each species sub-block.

The block builds, per species: the grouped array variable, its initial condition, a time
derivative kernel, and the self/opposite-absorption kernel
([GenericVacancyVoidGroupedNodalKernel.md] or [GenericInterstitialLoopGroupedNodalKernel.md]
depending on `geometry`). When both species are present and
[!param](/DefectDynamics/DefectDynamicsAction/opposite_absorption) is true (the default),
opposite-species absorption is active -- every "opposite" parameter the two kernels need (the
other species' diffusivity, dislocation bias, and grid) is derived from the other species' own
sub-block, never typed by the user. [PointDefectRecombinationNodalKernel.md] is added,
symmetrically, only when
[!param](/DefectDynamics/DefectDynamicsAction/recombination_radius) is given. The two are
independent switches: either can be on without the other.

The two species sub-blocks' `geometry` need not differ: `geometry = loop` for both is supported
(e.g. modeling a vacancy species as a loop instead of a void, following Pokor et al. 2004's own
Eqs 4a-4f, which use the same two-dimensional dislocation-bias absorption form for both species --
see [theory/defect_dynamics.md]). `geometry = void` for both species is not supported. When both
species are loops, the action also tells each loop kernel that its partner is a loop
(`opposite_geometry = loop`, with the partner's Burgers vector and bias), so that the drain on each
species' monomer by the other species' loops uses the loop absorption coefficient, the same one the
other kernel uses to shrink its loops. A void partner keeps the void coefficient.

The block also derives several quantities that would otherwise have to be computed outside MOOSE
and typed in by hand:

- Each species' dislocation sink coefficient, `sink = D*Z*dislocation_density`.
- The recombination rate, `recombination_rate = 4*pi*recombination_radius*(D_i + D_v)/atomic_volume`.
- Each species' capillary formation-energy coefficient and constant. By default,
  `A = (Ef - Eb2)/(2^exponent - 1)`, `C = Ef - A` from
  [!param](/DefectDynamics/DefectSpeciesAction/Ef) and
  [!param](/DefectDynamics/DefectSpeciesAction/Eb2); a species may instead give `A`/`C` directly
  via [!param](/DefectDynamics/DefectSpeciesAction/capillary_coefficient_eV) and
  [!param](/DefectDynamics/DefectSpeciesAction/capillary_constant_eV) -- see
  [DefectSpeciesAction.md] for when to use which, and why they are mutually exclusive.
- Each species' initial monomer concentration, `4*exp(-Ef/kT)` (using `Ef` from whichever of the
  two capillary input paths above was used).
- Each species' maximum cluster size (`num_cluster_sizes`), when given as
  [!param](/DefectDynamics/DefectSpeciesAction/max_cluster_radius_m) instead of directly.
- If a species sets [!param](/DefectDynamics/DefectSpeciesAction/resolvable_radius), the
  physical cluster size `n_minimum` and monomer radius `r1` needed by a number-density and
  average-radius postprocessor pair counting only clusters at or above that radius -- the
  smallest cluster size resolvable by whatever characterization technique (e.g. TEM) the input
  is being compared against. This is independent per species, not shared.

## Fokker-Planck Continuum Tail

Instead of grouping the large sizes in bins, the block can store them as a Fokker-Planck continuum
(Pokor et al. 2004, Section 4.1) by setting [!param](/DefectDynamics/DefectDynamicsAction/tail)
`= fokker_planck`. Sizes up to [!param](/DefectDynamics/DefectDynamicsAction/switch_size) are then
individual unknowns and larger sizes are finite-volume cells that start at `switch_size` + 1/2 with
the width [!param](/DefectDynamics/DefectDynamicsAction/continuum_initial_width) and grow by
[!param](/DefectDynamics/DefectDynamicsAction/continuum_width_ratio) from one cell to the next. Each
species builds a [ContinuumLoopNodalKernel.md] in place of the grouped loop kernel, the monomer and
the resolvable-population postprocessors read the continuum
([ClusterContinuumMoment.md]), and the size distribution is reported by
[ClusterContinuumDistribution.md]. Cascade seeding, `[Recovery]`, recombination, opposite-species
absorption and the constant dislocation sink are available as with grouping.
[!param](/DefectDynamics/DefectDynamicsAction/continuum_limiter) selects the second-order
(minmod-limited) or the first-order face values.

The continuum reaches a given accuracy with fewer unknowns than grouping, but it is narrower:

- Only the capillary binding-energy model is supported; a species cannot use a binding- or
  formation-energy table.
- Both species must have the same maximum cluster size, since each family reads the cells of the
  other. The grouping parameters of a species (`use_grouping`, `explicit_max`, `group_*`) are not
  allowed with the continuum, and `switch_size` must be at least 4 when a species has cascade
  seeding.
- The flux out of the last cell is zero, so the maximum cluster size must be large enough that the
  population does not reach it; check that the share of the point defects in the last cells is
  small.

In the Pokor benchmark at 603 K and 1 dpa with a largest size of 40,000, a switch size of 100 (576
unknowns for the two species) gives, at the end of the run, monomer concentrations, densities and
mean radii within 0.01 percent of the SUNDIALS Pokor benchmark (MNS_CD) with the same physics, in
about 10 seconds and 660 MB. A grouping with 968 unknowns and 1.7 GB comes within 10 percent of the
same reference, and groupings with 356 to 824 unknowns give resolvable radii that are too large by
factors of 2 to 6.
The time to reach the first resolvable loops is a steep onset, so the time step must be limited
(for example with a `timestep_limiting_postprocessor` of $0.01\,t + 0.5$ s) to resolve it.

## Time Step for the Onset of the Resolvable Population

The resolvable population (the clusters above the smallest radius a microscope can see,
[!param](/DefectDynamics/DefectSpeciesAction/resolvable_radius)) rises by many orders of magnitude within a small
interval of dose, from a density near zero to the reporting floor of $10^{12}$ m$^{-3}$ and beyond. A time step that is
fixed by the time, such as `IterationAdaptiveDT` with a cap of $0.01\,t + 0.5$ s, cannot follow it: in the six Yang, MD and Pokor
cases the density was off by up to a factor of 60 at the onset, whichever of grouping and `tail = fokker_planck` was used.

Setting [!param](/DefectDynamics/DefectDynamicsAction/timestep_relative_change) adds the postprocessor
`<block name>_timestep_limit` ([RelativeChangeTimestepLimit.md]), which limits the relative change of the resolvable densities of the species
with a `resolvable_radius` to that fraction per step (0.1 is a good start; smaller values follow the rise more closely and take more steps), bounded above by
[!param](/DefectDynamics/DefectDynamicsAction/timestep_cap_fraction) $\times\, t\,+$
[!param](/DefectDynamics/DefectDynamicsAction/timestep_cap_minimum). The input then names it in the time stepper:

```
[Executioner]
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-5
    optimal_iterations = 6
    iteration_window = 1
    growth_factor = 1.3
    timestep_limiting_postprocessor = DefectDynamics_timestep_limit
  []
[]
```

The action cannot set the time stepper itself. The name of the postprocessor is the block name followed by `_timestep_limit`, and the block needs
a `resolvable_radius` on at least one species.

## Dislocation Recovery

By default, each species' dislocation sink uses a constant
[!param](/DefectDynamics/DefectDynamicsAction/dislocation_density). Adding an optional
`[Recovery]` sub-block (`DefectRecoveryAction`) replaces this with a time-evolving dislocation
density instead, following the closed-form solution of Pokor et al. (2004)'s dislocation-climb
model (see [theory/defect_dynamics.md] for the derivation):

```
[Recovery]
  rho0 = 1e15             # initial dislocation density [m^-2]
  K = 5e13                # recovery rate constant [m^-1*s^-1]
  burgers_vector = 2.08e-10   # Burgers vector of the dislocation network [m]
[]
```

When `[Recovery]` is given, the block adds one [GenericDislocationRecoveryNodalKernel.md]
instance per species, each using that species' own `D0`/`Em`/`Z` but the SAME shared
`rho0`/`K`/`burgers_vector` -- Pokor's own equations use one un-subscripted dislocation density
for both species' sink terms, not an independent one per species.
[!param](/DefectDynamics/DefectDynamicsAction/dislocation_density) is left unset in this case
(it becomes unused, since `[Recovery]` supplies the entire sink contribution instead of adding
to it); setting it to a nonzero value alongside `[Recovery]` triggers a warning rather than
silently being ignored. `[Recovery]`'s own `burgers_vector` is a property of the dislocation
network itself, distinct from any loop species' own geometric
[!param](/DefectDynamics/DefectSpeciesAction/burgers_vector).

## Cascade Spectrum

By default, each species' entire generation rate goes to the monomer row. A species turns on
direct seeding of cluster sizes 2, 3, and 4 simply by setting a nonzero
[!param](/DefectDynamics/DefectSpeciesAction/cascade_f2),
[!param](/DefectDynamics/DefectSpeciesAction/cascade_f3), or
[!param](/DefectDynamics/DefectSpeciesAction/cascade_f4) on its own sub-block, following Pokor
et al. (2004) Eqs 1a-1e (see [theory/defect_dynamics.md] for the conservation-convention
derivation):

```
[vacancy]
  ...
  cascade_f2 = 0.2   # fraction of generation seeded directly as dimers
  cascade_f3 = 0.05  # fraction seeded as cluster size 3
  cascade_f4 = 0.02  # fraction seeded as cluster size 4
[]
```

For each species whose `cascade_f2`/`cascade_f3`/`cascade_f4` sum to more than zero, the block
adds one [GenericCascadeSpectrumNodalKernel.md] instance and reduces that species' own monomer
`generation` term by the same fractions so total point-defect production stays conserved. The
two species are independent -- one can seed cascade sizes while the other stays monomer-only.

## Input Options at a Glance

Several quantities in this block can be given in more than one way. Each row below is a
mutually exclusive choice -- exactly one option in the row must be given (except
`resolvable_radius`/`opposite_absorption`/`recombination_radius`, which are plain optional
switches, each defaulting to "off").

| What | Option A | Option B | Scope |
|---|---|---|---|
| Capillary formation-energy law | [!param](/DefectDynamics/DefectSpeciesAction/Ef) + [!param](/DefectDynamics/DefectSpeciesAction/Eb2) (default; derives $A$, $C$) | [!param](/DefectDynamics/DefectSpeciesAction/capillary_coefficient_eV) + [!param](/DefectDynamics/DefectSpeciesAction/capillary_constant_eV) ($A$, $C$ given directly) | per species |
| Maximum cluster size | [!param](/DefectDynamics/DefectSpeciesAction/num_cluster_sizes) (a point-defect count, given directly) | [!param](/DefectDynamics/DefectSpeciesAction/max_cluster_radius_m) (a physical radius in meters; converted internally) | per species |
| Dislocation bias (loop only) | [!param](/DefectDynamics/DefectSpeciesAction/bias) `= size_dependent` (default; uses `Z`, `Z_loop`, `alpha`) | `bias = constant` (uses `Z` only; the action sets `Z_loop = Z / sqrt(b/(8 pi a0))` so the bias equals `Z` at every size) | per species |
| Storage of the large cluster sizes | [!param](/DefectDynamics/DefectDynamicsAction/tail) `= grouping` (default; per species `explicit_max` and the `group_*` parameters) | `= fokker_planck` (`switch_size` and the `continuum_*` parameters) | shared |
| Emission binding energy | [!param](/DefectDynamics/DefectSpeciesAction/binding_energy_model) `= capillary` (default; always uses the law above) | `= binding_energy_table` or `= formation_energy_table` (direct DFT/MD data for small $n$, falling back to the capillary law beyond the table) | per species |

Independent on/off switches (not mutually exclusive with anything, each off by default):

| Switch | Effect when set | Scope |
|---|---|---|
| [!param](/DefectDynamics/DefectSpeciesAction/resolvable_radius) `> 0` | Adds a number-density/average-radius postprocessor pair counting only clusters at or above this radius | per species |
| [!param](/DefectDynamics/DefectDynamicsAction/recombination_radius) `> 0` | Adds the two symmetric point-defect recombination kernels (requires both species) | shared |
| [!param](/DefectDynamics/DefectDynamicsAction/opposite_absorption) `= false` | Turns off opposite-species absorption, independently of recombination (requires both species) | shared |
| [!param](/DefectDynamics/DefectSpeciesAction/cascade_f2)/[!param](/DefectDynamics/DefectSpeciesAction/cascade_f3)/[!param](/DefectDynamics/DefectSpeciesAction/cascade_f4), sum `> 0` | Seeds that species' cluster sizes 2-4 directly from the cascade, reducing its own monomer generation term correspondingly | per species |

"Per species" means the choice is made separately for `[vacancy]` and `[interstitial]` -- one
species can use `Ef`/`Eb2` while the other uses `capillary_coefficient_eV`/
`capillary_constant_eV`, for instance, or one can set `num_cluster_sizes` while the other sets
`max_cluster_radius_m`. "Shared" means the one top-level `[DefectDynamics]` block controls it
for the whole system.

## Direct Binding- or Formation-Energy Data

Each species sub-block also accepts
[!param](/DefectDynamics/DefectSpeciesAction/binding_energy_model),
[!param](/DefectDynamics/DefectSpeciesAction/binding_energy_table_eV), and
[!param](/DefectDynamics/DefectSpeciesAction/formation_energy_table_eV), forwarded unchanged to
the underlying kernel's own parameters of the same name. See
[GenericVacancyVoidGroupedNodalKernel.md] (void) or
[GenericInterstitialLoopGroupedNodalKernel.md] (loop) for how these select between the capillary
law and a direct DFT/MD-sourced binding- or formation-energy table for small cluster sizes, with
an automatic fallback to the capillary law beyond the table's range.

## Full Cluster Size Distribution

By default, only aggregate scalars are available per species
([!param](/DefectDynamics/DefectSpeciesAction/resolvable_radius)'s postprocessor pair, and the
always-added monomer concentration). A species additionally sets
[!param](/DefectDynamics/DefectSpeciesAction/output_size_distribution) `= true` to get its full
per-size concentration (and radius) distribution -- one row per physical cluster size $n$ from 1
to that species' own `num_cluster_sizes`, reconstructed via [ClusterSizeDistributionGrouped.md]
the same way the grouped kernels reconstruct a cluster size's own concentration internally:

```
[DefectDynamics]
  ...
  size_distribution_times = '1.7e7 3.4e7'   # simulation times [s], NOT dose/dpa, shared by every species

  [vacancy]
    ...
    output_size_distribution = true
  []
[]
```

[!param](/DefectDynamics/DefectDynamicsAction/size_distribution_times) is shared by every
species (not per-species) so that species with this output turned on stay directly comparable at
the same instants, and must be non-empty if any species turns the output on. The block adds a
dedicated `CSV` output, restricted to exactly these times (`sync_times`/`sync_only`) and to
vector postprocessors only, so the regular per-timestep scalar postprocessor history from the
general `[Outputs]` block is completely unaffected; MOOSE's own naming convention produces one
file per species per requested time, e.g.
`<file_base>_<DefectDynamics-block-name>_size_distribution_csv_vacancy_size_distribution_0002.csv`
for a requested time reached at time step 2 (the number in the name is the time step, not the
position of the time in the list). Within a grouped (large-$n$) bin, the reconstructed values lie on that
bin's own linear fit (its stored $L_0$/$L_1$ moments), so they are not independent data -- see
[ClusterSizeDistributionGrouped.md] for the reconstruction formula.

## Console Output at Large Cluster Counts

MOOSE's console prints a residual norm for the variables with the largest residuals, and it
computes that norm separately for every component of an array variable. With hundreds of
thousands of cluster sizes, that costs a noticeable share of each time step (about a quarter of the
run time in a case with 400,000 sizes per species). Set
`outlier_variable_norms = false` in the console output to skip it; the solution is not affected.

```
[Outputs]
  [console]
    type = Console
    outlier_variable_norms = false
  []
[]
```

## Disabling or Omitting a Term

Setting `dislocation_density = 0` turns off the dislocation sink for every species (unless
`[Recovery]` is given, in which case `dislocation_density` is unused and `[Recovery]`'s own
`rho0` controls the sink instead -- omitting the `[Recovery]` sub-block entirely is what turns
off the time-evolving sink). Omitting `recombination_radius` (or leaving it at its default of
zero) omits the recombination kernels.
Setting `opposite_absorption = false` turns off opposite-species absorption while leaving
recombination independently controlled by `recombination_radius` -- useful for isolating one
coupling physics from the other, or for running two species with no coupling between them at
all. Supplying only one species sub-block (`[vacancy]` or `[interstitial]` alone) reduces that
species' kernel to self-reaction-only behavior, with no opposite-species or recombination terms
at all -- the same self-reaction-only mode the low-level kernels already support when
`opposite_coupled_var` is left unset. Leaving a species'
[!param](/DefectDynamics/DefectSpeciesAction/resolvable_radius) at its default of zero omits
that species' number-density/average-radius postprocessor pair. Leaving a species' own
`cascade_f2`/`cascade_f3`/`cascade_f4` at their default of zero (or setting them back to zero)
is what turns off cascade seeding for that species.

## Example Input Syntax

!listing test/tests/defect_dynamics/point_defect_recombination_small_action.i block=DefectDynamics

!syntax parameters /DefectDynamics/DefectDynamicsAction

## See Also

- [GenericVacancyVoidGroupedNodalKernel.md]
- [GenericInterstitialLoopGroupedNodalKernel.md]
- [PointDefectRecombinationNodalKernel.md]
- [GenericDislocationRecoveryNodalKernel.md]
- [GenericCascadeSpectrumNodalKernel.md]
- [ClusterSizeDistributionGrouped.md]
- [ContinuumLoopNodalKernel.md]
- [ClusterContinuumMoment.md]
- [ClusterContinuumDistribution.md]
- [theory/defect_dynamics.md]
