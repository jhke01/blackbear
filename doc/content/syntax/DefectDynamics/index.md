# DefectDynamics

!syntax description /DefectDynamics/DefectDynamicsAction

## Description

The `DefectDynamics` block is a species-organized front end to the grouped defect-dynamics
kernels ([GenericVacancyVoidGroupedNodalKernel.md], [GenericInterstitialLoopGroupedNodalKernel.md],
[PointDefectRecombinationNodalKernel.md]; see [theory/defect_dynamics.md] for the underlying
equations). One or two named species sub-blocks (commonly `[vacancy]` and `[interstitial]`,
though the name is only a label) describe each species' own diffusivity, capillary
formation-energy law, and (for a loop) dislocation bias; one `[Grid]` sub-block describes the
grouped-tail layout shared by both species.

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
- If [!param](/DefectDynamics/DefectDynamicsAction/tem_visible_radius) is given, the physical
  cluster size `n_minimum` and monomer radius `r1` needed by a TEM-visible number-density and
  average-radius postprocessor pair, for both species.

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

Setting `dislocation_density = 0` turns off the dislocation sink for every species. Omitting
`recombination_radius` (or leaving it at its default of zero) omits the recombination kernels.
Setting `opposite_absorption = false` turns off opposite-species absorption while leaving
recombination independently controlled by `recombination_radius` -- useful for isolating one
coupling physics from the other, or for running two species with no coupling between them at
all. Supplying only one species sub-block (`[vacancy]` or `[interstitial]` alone) reduces that
species' kernel to self-reaction-only behavior, with no opposite-species or recombination terms
at all -- the same self-reaction-only mode the low-level kernels already support when
`opposite_coupled_var` is left unset.

## Example Input Syntax

!listing test/tests/defect_dynamics/point_defect_recombination_small_action.i block=DefectDynamics

!syntax parameters /DefectDynamics/DefectDynamicsAction

## See Also

- [GenericVacancyVoidGroupedNodalKernel.md]
- [GenericInterstitialLoopGroupedNodalKernel.md]
- [PointDefectRecombinationNodalKernel.md]
- [theory/defect_dynamics.md]
