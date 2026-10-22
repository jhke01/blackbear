# ADInterstitialLoopGroupedNodalKernel / InterstitialLoopGroupedNodalKernel

!syntax description /NodalKernels/InterstitialLoopGroupedNodalKernel

## Description

See [theory/defect_dynamics.md] for a full derivation connecting this
kernel's equations to [GenericVacancyVoidGroupedNodalKernel.md] and
[PointDefectRecombinationNodalKernel.md].

`InterstitialLoopGroupedNodalKernel` solves the self-reaction (absorption and
thermal emission) terms of an interstitial-loop population balance, using the
same grouped L0/L1 representation as
[GenericClusterDynamicsGroupedNodalKernel.md]: cluster sizes through
[!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/explicit_max) are
stored explicitly, and the remaining physical sizes are divided into bins
carrying two moments each. See that page for the grouped representation,
moment evolution equations, reconstruction limiter, and intra-variable
Jacobian, which this kernel reuses unchanged.

This kernel implements the interstitial family's own growth/shrinkage by
interstitial-monomer absorption and emission, plus the opposite-species
absorption terms described below: a vacancy monomer absorbed onto a loop
shrinks it by one interstitial. Point-defect recombination with a vacancy
variable and dislocation sink loss for the monomer (via the
[!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/sink) parameter) are
added separately by [PointDefectRecombinationNodalKernel.md] and the `sink`
parameter, respectively; mobile interstitial-cluster reactions are not yet
implemented.

## Rate Coefficients

Interstitial loops are two-dimensional clusters, so the loop radius follows

!equation
r_n = \sqrt{\frac{nV_{at}}{\pi b}},

where $b$ is the Burgers vector
([!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/burgers_vector),
defaulting to
[!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/lattice_parameter)$/\sqrt{3}$
when left at zero). This differs from the $n^{1/3}$ radius law used for
three-dimensional clusters in
[GenericClusterDynamicsNodalKernel.md] and
[GenericVacancyVoidGroupedNodalKernel.md].

The interstitial-monomer diffusivity follows an Arrhenius law,
$D_i = D_{0,i}\exp(-E_{m,i}/k_BT)$, and the absorption rate coefficient uses a
size-dependent dislocation-loop bias factor in place of the void kernel's
finite-sink-strength correction:

!equation
\beta_n = 2\pi r_n Z_i(n) D_i / V_{at},
\qquad
Z_i(n) = Z_{i,0} + \left[\sqrt{\frac{b}{8\pi a_0}}Z_{loop,i} - Z_{i,0}\right]
n^{-\alpha_{loop,i}/2},

where $a_0$ is
[!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/lattice_parameter),
$Z_{i,0}$, $Z_{loop,i}$, and $\alpha_{loop,i}$ are
[!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/z_i),
[!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/z_loop_i), and
[!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/alpha_loop_i). $Z_i(n)$
is $\sqrt{b/(8\pi a_0)}\,Z_{loop,i}$ at $n=1$ (about $0.152\,Z_{loop,i}$ for
$b=a_0/\sqrt3$) and approaches $Z_{i,0}$ as $n$ grows. So $Z_{i,0}$ is the large-loop
bias and $Z_{loop,i}$ sets the small-loop bias. A constant bias equal to
$Z_{i,0}$ at every size needs $Z_{loop,i}=Z_{i,0}/\sqrt{b/(8\pi a_0)}$.

Thermal emission follows detailed balance from a binding energy $E_b(n)$,

!equation
\alpha_n = \beta_{n-1}\exp\left(-\frac{E_b(n)}{k_BT}\right).

By default (and always as the fallback described below), $E_b(n)$ comes from
a capillary-law formation energy, with the same structure as
[GenericVacancyVoidGroupedNodalKernel.md]'s void emission model but the
loop's own exponent and coefficients:

!equation
E_f(n) = A_i n^{p_i} + C_i,
\qquad
E_b(n) = E_{form,i} - \left[E_f(n)-E_f(n-1)\right],

where $E_{form,i}$ is
[!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/formation_energy_eV).
These, together with $A_i$/$C_i$/$p_i$, are always required, in every
[!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/binding_energy_model)
setting -- see below.

## Direct Binding- or Formation-Energy Data

[!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/binding_energy_model)
selects how $E_b(n)$ is obtained, with the same three modes, table indexing,
and capillary-fallback behavior as
[GenericVacancyVoidGroupedNodalKernel.md]'s own
`binding_energy_model`/`binding_energy_table_eV`/`formation_energy_table_eV` --
see that page for the full discussion of why the fallback (rather than
[GenericClusterDynamicsNodalKernel.md]'s hard error past the table's end) is
needed at this model's scale.

Both table modes are exercised in

!listing test/tests/defect_dynamics/interstitial_loop_energy_models.i

using the real MD binding energies from Soneda and Diaz de la Rubia (1998),
*Philosophical Magazine A* 78:5, 995-1019 -- their Table 1, interstitial
clusters $n=2..10$ -- with
[!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/num_cluster_sizes) = 30
to also exercise the capillary fallback beyond $n=10$.

## Opposite-Species Absorption

When [!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/opposite_coupled_var)
is provided (a [GenericVacancyVoidGroupedNodalKernel.md] array variable), this
kernel also couples in the opposite species' effect on the loop population,
following `oppositeShrinkFluxAt` in the reference AMSS rate-theory model. For
a loop of size $n\geq 2$, an absorbed vacancy monomer shrinks it to size
$n-1$ at rate

!equation
\beta^{opp}_n = 2\pi r_n Z_v(n) D_v / V_{at},
\qquad
Z_v(n) = Z_{v,0} + \left[\sqrt{\frac{b}{8\pi a_0}}Z_{loop,v} - Z_{v,0}\right]
n^{-\alpha_{loop,v}/2},

the same functional form as $\beta_n$ above but with the vacancy-absorption
bias parameters
[!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/z_v),
[!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/z_loop_v), and
[!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/alpha_loop_v) in
place of the self-absorption bias, and the vacancy diffusivity
$D_v=D_{0,v}\exp(-E_{m,v}/k_BT)$
([!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/D0_v) and
[!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/Em_v_eV)) in place
of $D_i$. This shrinkage flux is subtracted from the same flux array that
[GenericClusterDynamicsGroupedNodalKernel.md]'s moment equations already use,
so no other structural change is needed for explicit sizes or grouped bins;
the monomer equation additionally gains a $+\beta^{opp}_2 C_{v,1} C_2$ term (one
unit returned to the monomer pool, not doubled the way self-dimerization is).

The corresponding loss from the vacancy monomer pool is computed by summing
the *entire* coupled array. Each term uses the coefficient at which the coupled family's own
clusters absorb this kernel's monomer, with this kernel's own $D_i$ (no new diffusivity parameter
is needed, since this kernel already computes $D_i$ for its self-absorption terms). Which
coefficient that is depends on the coupled family's geometry,
[!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/opposite_geometry):

- `void` (the default): the void family's finite-sink-strength coefficient, the same form as
  [GenericVacancyVoidGroupedNodalKernel.md]'s own $\beta_n$ but with $D_i$.
- `loop`: the loop coefficient $2\pi r_n Z(n) D_i / V_{at}$, with the coupled loops' radius law
  (their Burgers vector,
  [!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/opposite_burgers_vector), defaulting to
  this kernel's own) and their bias for capturing this kernel's monomer,
  [!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/opposite_z_i),
  [!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/opposite_z_loop_i), and
  [!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/opposite_alpha_loop_i). Use this when
  both species are loops, e.g. a vacancy modeled as a loop following Pokor et al. (2004).

This loss and the shrinkage of the coupled family's clusters are the two sides of one absorption
event, so they must use one rate. With the wrong geometry the monomer pool is drained at a
different rate from the one at which the clusters shrink, and defects are no longer conserved.
For two loop families with no generation, sinks or recombination, the total vacancy content minus
the total interstitial content stays constant; see
!listing test/tests/defect_dynamics/defect_dynamics_two_loop_conservation.i
and its Jacobian check
!listing test/tests/defect_dynamics/defect_dynamics_two_loop_jacobian.i
Reconstructing the coupled void array requires its own grouping layout,
described independently of this kernel's own layout by the
`opposite_num_cluster_sizes`/`opposite_use_grouping`/`opposite_explicit_max`/
`opposite_group_spacing`/`opposite_group_linear_width`/
`opposite_group_geometric_ratio`/`opposite_enforce_group_nonnegative`/
`opposite_group_nonnegative_tolerance_factor` parameters, which must match the
coupled kernel instance's own grouping
settings.

When `opposite_coupled_var` is not provided, every opposite-absorption term
evaluates to exactly zero and this kernel reduces to the self-reaction-only
behavior described above.

!alert note
The opposite-absorption terms are implemented only for the non-AD
`InterstitialLoopGroupedNodalKernel` variant, for the same reason as
[PointDefectRecombinationNodalKernel.md]: `Coupleable::coupledArrayValue()`
returns a plain, non-AD array value, and `ADArrayNodalKernel::setJacobian()`
unconditionally errors, so the AD variant's cross-variable Jacobian block
would silently assemble to zero. This kernel's self-reaction terms remain
AD-correct in `ADInterstitialLoopGroupedNodalKernel`; only the
opposite-absorption Jacobian contributions require the non-AD variant.

## Required and Recommended Solver Settings

This kernel requires `use_hash_table_matrix_assembly = true` in the
`[Problem]` block, for the same reason as
[GenericClusterDynamicsGroupedNodalKernel.md]. See that page for the full
discussion of recommended `[Preconditioning]` and `[Executioner]` settings.

!alert note
Realistic Arrhenius diffusivities make $\beta_1$ many orders of magnitude
larger than a typical PETSc finite-difference Jacobian test can resolve with
one global step length (the `wp` algorithm used by `-snes_test_jacobian`).
`test/tests/defect_dynamics/interstitial_loop_jacobian_grouped.i` uses a
strongly reduced test-only `D0_i` for this reason; see the comments in that
file.

## Example Input Syntax

!listing test/tests/defect_dynamics/interstitial_loop_grouped_small.i

The opposite-species absorption terms are exercised in
!listing test/tests/defect_dynamics/defect_dynamics_opposite_absorption_small.i

!syntax parameters /NodalKernels/InterstitialLoopGroupedNodalKernel

!syntax inputs /NodalKernels/InterstitialLoopGroupedNodalKernel

!syntax children /NodalKernels/InterstitialLoopGroupedNodalKernel

## See Also

- [theory/defect_dynamics.md]
- [GenericClusterDynamicsGroupedNodalKernel.md]
- [GenericVacancyVoidGroupedNodalKernel.md]
- [PointDefectRecombinationNodalKernel.md]
- [GroupedArrayFunctionIndexIC.md]
