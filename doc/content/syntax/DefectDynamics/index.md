# DefectDynamics

!syntax description /DefectDynamics/DefectDynamicsAction

## Description

The `DefectDynamics` block is a species-organized front end to the grouped defect-dynamics
kernels ([GenericVacancyVoidGroupedNodalKernel.md], [GenericInterstitialLoopGroupedNodalKernel.md],
[PointDefectRecombinationNodalKernel.md], [GenericDislocationRecoveryNodalKernel.md],
[GenericCascadeSpectrumNodalKernel.md]; see [theory/defect_dynamics.md] for the underlying
equations). One or two named species sub-blocks (commonly `[vacancy]` and `[interstitial]`,
though the name is only a label) describe each species' own diffusivity, capillary
formation-energy law, (for a loop) dislocation bias, and its own maximum cluster size; one
`[Grid]` sub-block describes the grouped-tail *shape* (explicit size, bin spacing, growth ratio)
shared by both species, and one optional `[Recovery]` sub-block describes a dislocation network
shared by both species whose density evolves in time instead of staying fixed. A species also
turns on cascade seeding of cluster sizes 2-4 simply by setting a nonzero cascade fraction on its
own sub-block.

Void and loop populations routinely grow to very different physical sizes under the same
irradiation conditions, so the maximum cluster size itself is deliberately NOT part of the
shared `[Grid]` block -- each species sets its own via
[!param](/DefectDynamics/DefectSpeciesAction/num_cluster_sizes) (the cap given directly, as a
point-defect count) or
[!param](/DefectDynamics/DefectSpeciesAction/max_cluster_radius_m) (the cap given as a physical
radius in meters, converted internally using that species' own monomer radius law); exactly one
of the two must be given per species.

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
| Dislocation bias (loop only) | [!param](/DefectDynamics/DefectSpeciesAction/bias) `= size_dependent` (default; uses `Z`, `Z_loop`, `alpha`) | `bias = constant` (uses `Z` only, equivalent to `Z_loop = Z`) | per species |
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
- [theory/defect_dynamics.md]
