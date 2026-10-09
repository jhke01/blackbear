# ADInterstitialLoopGroupedNodalKernel / InterstitialLoopGroupedNodalKernel

!syntax description /NodalKernels/InterstitialLoopGroupedNodalKernel

## Description

`InterstitialLoopGroupedNodalKernel` solves the self-reaction (absorption and
thermal emission) terms of an interstitial-loop population balance, using the
same grouped L0/L1 representation as
[GenericClusterDynamicsGroupedNodalKernel.md]: cluster sizes through
[!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/explicit_max) are
stored explicitly, and the remaining physical sizes are divided into bins
carrying two moments each. See that page for the grouped representation,
moment evolution equations, reconstruction limiter, and intra-variable
Jacobian, which this kernel reuses unchanged.

This kernel implements only the interstitial family's own growth/shrinkage by
interstitial-monomer absorption and emission. Point-defect recombination with
a vacancy variable and dislocation sink loss for the monomer (via the
[!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/sink) parameter) are
the only inter-family or sink physics present; mobile interstitial-cluster
reactions are not yet implemented.

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
interpolates from the small-loop bias $Z_{i,0}$ to the large-loop bias
$Z_{loop,i}$ as $n$ grows.

Thermal emission follows detailed balance from a capillary-law formation
energy, with the same structure as
[GenericVacancyVoidGroupedNodalKernel.md]'s void emission model but the
loop's own exponent and coefficients:

!equation
E_f(n) = A_i n^{p_i} + C_i,
\qquad
E_b(n) = E_{form,i} - \left[E_f(n)-E_f(n-1)\right],
\qquad
\alpha_n = \beta_{n-1}\exp\left(-\frac{E_b(n)}{k_BT}\right).

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

!syntax parameters /NodalKernels/InterstitialLoopGroupedNodalKernel

!syntax inputs /NodalKernels/InterstitialLoopGroupedNodalKernel

!syntax children /NodalKernels/InterstitialLoopGroupedNodalKernel

## See Also

- [GenericClusterDynamicsGroupedNodalKernel.md]
- [GenericVacancyVoidGroupedNodalKernel.md]
- [GroupedArrayFunctionIndexIC.md]
