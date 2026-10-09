# ADVacancyVoidGroupedNodalKernel / VacancyVoidGroupedNodalKernel

!syntax description /NodalKernels/VacancyVoidGroupedNodalKernel

## Description

`VacancyVoidGroupedNodalKernel` solves the self-reaction (absorption and thermal
emission) terms of a vacancy-void population balance, using the same grouped
L0/L1 representation as [GenericClusterDynamicsGroupedNodalKernel.md]: cluster
sizes through
[!param](/NodalKernels/VacancyVoidGroupedNodalKernel/explicit_max) are stored
explicitly, and the remaining physical sizes are divided into bins carrying two
moments each. See that page for the grouped representation, moment evolution
equations, reconstruction limiter, and intra-variable Jacobian, which this
kernel reuses unchanged.

This kernel implements the vacancy family's own growth/shrinkage by
vacancy-monomer absorption and emission, plus the opposite-species absorption
terms described below: an interstitial monomer absorbed onto a void shrinks
it by one vacancy. Point-defect recombination with an interstitial variable
and dislocation sink loss for the monomer (via the
[!param](/NodalKernels/VacancyVoidGroupedNodalKernel/sink) parameter) are
added separately by [PointDefectRecombinationNodalKernel.md] and the `sink`
parameter, respectively; mobile interstitial-cluster reactions are not yet
implemented.

## Rate Coefficients

Voids are three-dimensional clusters, so the void radius follows

!equation
r_n = \left(\frac{3nV_{at}}{4\pi}\right)^{1/3},

the same functional form as the spherical-cluster radius in
[GenericClusterDynamicsNodalKernel.md]. The vacancy-monomer diffusivity follows
an Arrhenius law,

!equation
D_v = D_{0,v}\exp\left(-\frac{E_{m,v}}{k_BT}\right),

and the absorption rate coefficient includes a finite-sink-strength
correction that does not appear in the Cu-precipitation kernel:

!equation
\beta_n = \frac{4\pi r_n D_v}{\left(1+\dfrac{a_0}{r_n}\right)V_{at}},

where $a_0$ is
[!param](/NodalKernels/VacancyVoidGroupedNodalKernel/lattice_parameter). This
correction becomes significant once $r_n$ approaches $a_0$ and suppresses the
absorption rate relative to the uncorrected $4\pi r_n D_v/V_{at}$ form.

Thermal emission follows detailed balance from a capillary-law formation
energy,

!equation
E_f(n) = A_v n^{p_v} + C_v,
\qquad
E_b(n) = E_{form,v} - \left[E_f(n)-E_f(n-1)\right],
\qquad
\alpha_n = \beta_{n-1}\exp\left(-\frac{E_b(n)}{k_BT}\right),

where $p_v$ is
[!param](/NodalKernels/VacancyVoidGroupedNodalKernel/void_energy_exponent),
$A_v$ and $C_v$ are
[!param](/NodalKernels/VacancyVoidGroupedNodalKernel/void_energy_coefficient_eV)
and
[!param](/NodalKernels/VacancyVoidGroupedNodalKernel/void_energy_constant_eV),
and $E_{form,v}$ is
[!param](/NodalKernels/VacancyVoidGroupedNodalKernel/formation_energy_eV).

## Opposite-Species Absorption

When [!param](/NodalKernels/VacancyVoidGroupedNodalKernel/opposite_coupled_var)
is provided (an [GenericInterstitialLoopGroupedNodalKernel.md] array variable), this
kernel also couples in the opposite species' effect on the vacancy population,
following `oppositeShrinkFluxAt` in the reference AMSS rate-theory model. For
a void of size $n\geq 2$, an absorbed interstitial monomer shrinks it to size
$n-1$ at rate

!equation
\beta^{opp}_n = \frac{4\pi r_n D_i}{\left(1+\dfrac{a_0}{r_n}\right)V_{at}},

the same finite-sink-strength form as $\beta_n$ above, but with the
interstitial diffusivity $D_i=D_{0,i}\exp(-E_{m,i}/k_BT)$ in place of $D_v$
([!param](/NodalKernels/VacancyVoidGroupedNodalKernel/D0_i) and
[!param](/NodalKernels/VacancyVoidGroupedNodalKernel/Em_i_eV)). This shrinkage
flux is subtracted from the same flux array that
[GenericClusterDynamicsGroupedNodalKernel.md]'s moment equations already use,
so no other structural change is needed for explicit sizes or grouped bins;
the monomer equation additionally gains a $+\beta^{opp}_2 C_{i,1} C_2$ term (one
unit returned to the monomer pool, not doubled the way self-dimerization is).

The corresponding loss from the interstitial monomer pool is computed by
summing the *entire* coupled loop array using the loop family's own
absorption-of-vacancy coefficient,

!equation
\beta^{opp}_{loop,n} = 2\pi r_{loop,n} Z_v(n) D_v / V_{at},

where $Z_v(n)$ has the same functional form as
[GenericInterstitialLoopGroupedNodalKernel.md]'s own bias factor but with the
vacancy-side parameters
[!param](/NodalKernels/VacancyVoidGroupedNodalKernel/opposite_z_v),
[!param](/NodalKernels/VacancyVoidGroupedNodalKernel/opposite_z_loop_v), and
[!param](/NodalKernels/VacancyVoidGroupedNodalKernel/opposite_alpha_loop_v).
Reconstructing the coupled loop array requires its own grouping layout,
described independently of this kernel's own layout by the
`opposite_num_cluster_sizes`/`opposite_use_grouping`/`opposite_explicit_max`/
`opposite_group_spacing`/`opposite_group_linear_width`/
`opposite_group_geometric_ratio`/`opposite_enforce_group_nonnegative`/
`opposite_group_nonnegative_tolerance_factor` parameters, which must match the
coupled [GenericInterstitialLoopGroupedNodalKernel.md] instance's own grouping
settings.

When `opposite_coupled_var` is not provided, every opposite-absorption term
evaluates to exactly zero and this kernel reduces to the self-reaction-only
behavior described above.

!alert note
The opposite-absorption terms are implemented only for the non-AD
`VacancyVoidGroupedNodalKernel` variant, for the same reason as
[PointDefectRecombinationNodalKernel.md]: `Coupleable::coupledArrayValue()`
returns a plain, non-AD array value, and `ADArrayNodalKernel::setJacobian()`
unconditionally errors, so the AD variant's cross-variable Jacobian block
would silently assemble to zero. This kernel's self-reaction terms remain
AD-correct in `ADVacancyVoidGroupedNodalKernel`; only the opposite-absorption
Jacobian contributions require the non-AD variant.

## Required and Recommended Solver Settings

This kernel requires `use_hash_table_matrix_assembly = true` in the
`[Problem]` block, for the same reason as
[GenericClusterDynamicsGroupedNodalKernel.md]. See that page for the full
discussion of recommended `[Preconditioning]` and `[Executioner]` settings.

!alert note
Realistic Arrhenius diffusivities make $\beta_1$/$\beta^{opp}_1$ many orders of
magnitude larger than a typical PETSc finite-difference Jacobian test can
resolve with one global step length (the `wp` algorithm used by
`-snes_test_jacobian`).
`test/tests/defect_dynamics/defect_dynamics_opposite_absorption_jacobian.i`
uses strongly reduced test-only `D0_v`/`D0_i` for this reason; see the
comments in that file.

## Example Input Syntax

!listing test/tests/defect_dynamics/vacancy_void_grouped_small.i

The opposite-species absorption terms are exercised in
!listing test/tests/defect_dynamics/defect_dynamics_opposite_absorption_small.i

!syntax parameters /NodalKernels/VacancyVoidGroupedNodalKernel

!syntax inputs /NodalKernels/VacancyVoidGroupedNodalKernel

!syntax children /NodalKernels/VacancyVoidGroupedNodalKernel

## See Also

- [GenericClusterDynamicsGroupedNodalKernel.md]
- [GenericInterstitialLoopGroupedNodalKernel.md]
- [PointDefectRecombinationNodalKernel.md]
- [GroupedArrayFunctionIndexIC.md]
