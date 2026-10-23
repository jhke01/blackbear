# ContinuumLoopNodalKernel

!syntax description /NodalKernels/ContinuumLoopNodalKernel

## Description

`ContinuumLoopNodalKernel` evolves the size distribution of one family of two-dimensional loops or
three-dimensional voids (set by [!param](/NodalKernels/ContinuumLoopNodalKernel/geometry)) with a
representation that is discrete for small sizes and a Fokker-Planck continuum for large ones, as
described for loops larger than about 100 atoms by [!cite](pokor2004). It is an alternative to the
grouped $L_0/L_1$ representation of [GenericInterstitialLoopGroupedNodalKernel.md] and
[GenericVacancyVoidGroupedNodalKernel.md]: the physics (the radius, the size-dependent bias of a
loop, the finite-sink-strength absorption of a void, detailed-balance emission from a capillary-law
binding energy) is the same, and only the way the large sizes are stored differs. A void ignores
[!param](/NodalKernels/ContinuumLoopNodalKernel/Z_loop) and
[!param](/NodalKernels/ContinuumLoopNodalKernel/alpha); its
[!param](/NodalKernels/ContinuumLoopNodalKernel/Z) only sets the large-loop limit of the bias with
which loops of the opposite family absorb its monomer
([!param](/NodalKernels/ContinuumLoopNodalKernel/opposite_geometry) selects the form of the
coefficient for the opposite family).

The kernel can model one loop family with a monomer source and a linear sink, or two families that
shrink by absorbing each other's monomers (see [Two coupled families](#two-families)). Recombination,
the cascade source and a time-dependent sink are separate kernels
([PointDefectRecombinationNodalKernel.md], [GenericCascadeSpectrumNodalKernel.md],
[GenericDislocationRecoveryNodalKernel.md]), as for the grouped kernels. The
[DefectDynamics](syntax/DefectDynamics/index.md) block builds all of them with
[!param](/DefectDynamics/DefectDynamicsAction/tail) `= fokker_planck`.

## Representation

Sizes $n = 1, \ldots, n_s$ ([!param](/NodalKernels/ContinuumLoopNodalKernel/switch_size)) are
individual unknowns. Above $n_s$ the concentration per unit size is averaged over finite-volume
cells that start at $n_s + 1/2$ with the width
[!param](/NodalKernels/ContinuumLoopNodalKernel/continuum_initial_width), grow by the factor
[!param](/NodalKernels/ContinuumLoopNodalKernel/continuum_width_ratio) from one cell to the next,
and end at $N + 1/2$, where $N$ is [!param](/NodalKernels/ContinuumLoopNodalKernel/num_cluster_sizes).
The array variable therefore has $n_s$ components for the monomer and the sizes up to $n_s$, plus one
component for each cell. A switch size of 100 and $N = 4 \times 10^6$ give 406 components.

Between two discrete sizes the flux is $J_n = \beta_n C_1 C_n - \alpha_{n+1} C_{n+1}$. Between two
cells it is the sum of an advective part, the drift $\beta C_1 - \alpha$ times the upwind face value,
and a diffusive part that acts on the cell averages with the coefficient $(\beta C_1 + \alpha)/2$.
The coefficients at a real-valued size are interpolated linearly between the neighboring integer
sizes. The flux out of the last cell is zero.

With [!param](/NodalKernels/ContinuumLoopNodalKernel/use_limiter) true the face values are
reconstructed from the cell averages and a limited slope (the generalized minmod limiter with
$\theta = 1.5$), which is second order and is the scheme of the SUNDIALS Pokor benchmark. With
`use_limiter = false` the face values are the cell averages, which is first-order upwind.

The Jacobian is written out by hand. With the limiter on, it includes the dependence of the slope on
the two neighboring cells, taking the derivative of the one-sided estimate that minmod selects. It is
exact except at the points where the selection changes.

## Two coupled families id=two-families

With [!param](/NodalKernels/ContinuumLoopNodalKernel/opposite_coupled_var), the loops also shrink by
absorbing a monomer of the other point-defect type, and the monomer balance of each family loses
what the other family's loops absorb (Pokor et al. 2004, Eqs. 2 and 3). The shrinking flux
$-\beta^{\mathrm{opp}}_{n+1} C^{\mathrm{opp}}_1 C_{n+1}$ has its own upwinding, since its drift is
always toward smaller sizes, and shares the limited slopes of the self-reaction flux. The two
families must have the same layout (the same
[!param](/NodalKernels/ContinuumLoopNodalKernel/num_cluster_sizes),
[!param](/NodalKernels/ContinuumLoopNodalKernel/switch_size) and cell widths).

The absorption coefficient of this family's loops for the opposite monomer is $2 \pi r_n Z_{\mathrm{opp}}(n)
D_{\mathrm{opp}} / V_{at}$, with the radius of this family's loops and the diffusivity of the
opposite monomer, set by [!param](/NodalKernels/ContinuumLoopNodalKernel/opposite_D0),
[!param](/NodalKernels/ContinuumLoopNodalKernel/opposite_Em_eV),
[!param](/NodalKernels/ContinuumLoopNodalKernel/opposite_Z),
[!param](/NodalKernels/ContinuumLoopNodalKernel/opposite_Z_loop) and
[!param](/NodalKernels/ContinuumLoopNodalKernel/opposite_alpha). The opposite family's loops absorb
this family's monomer with the radius set by
[!param](/NodalKernels/ContinuumLoopNodalKernel/opposite_burgers_vector) and a bias that approaches
[!param](/NodalKernels/ContinuumLoopNodalKernel/Z) with
[!param](/NodalKernels/ContinuumLoopNodalKernel/partner_Z_loop) and
[!param](/NodalKernels/ContinuumLoopNodalKernel/partner_alpha) (default: this family's own
`Z_loop` and `alpha`). The two kernels must use one set of values for the same absorption event,
which the [DefectDynamics](syntax/DefectDynamics/index.md) block ensures.

The Jacobian between the two variables is written out too, so the preconditioner needs
`off_diag_row` and `off_diag_column` for the two variables, and the memory warning for coupled array
variables applies.

## Choosing the largest size and the time step

The flux out of the last cell is zero, so a population that reaches the last cell piles up there.
In a test with $N = 40{,}000$, 99 percent of the point defects were in the last five cells after
$10^6$ s; with $N = 4 \times 10^6$ the fraction was $3 \times 10^{-9}$. Check the fraction of the
population in the last cells before relying on a result.

With `IterationAdaptiveDT` alone the step can grow faster than the transient can be followed, and the
monomer can end up off by tens of percent while every step converges. Limiting the step in proportion
to the time (a `timestep_limiting_postprocessor` equal to $0.01\,t + 0.5$ s in the tests) removes this.

## Example Input Syntax

In this example the loop family is discrete up to size 20 and a continuum up to size 4000, starting
from an empty distribution with a constant monomer source and a linear sink.

!listing test/tests/continuum_tail/continuum_loop_transient.i block=NodalKernels

In this example two coupled families, each discrete up to size 20 and a continuum up to size 2000,
shrink by absorbing each other's monomers. The kernel of the vacancy loops is shown.

!listing test/tests/continuum_tail/continuum_pair_jacobian.i block=NodalKernels/vacancy_dyn

!syntax parameters /NodalKernels/ContinuumLoopNodalKernel

!syntax inputs /NodalKernels/ContinuumLoopNodalKernel

!syntax children /NodalKernels/ContinuumLoopNodalKernel
