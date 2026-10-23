# ContinuumLoopNodalKernel

!syntax description /NodalKernels/ContinuumLoopNodalKernel

## Description

`ContinuumLoopNodalKernel` evolves the size distribution of one family of two-dimensional loops
with a representation that is discrete for small sizes and a Fokker-Planck continuum for large
ones, as described for loops larger than about 100 atoms by [!cite](pokor2004). It is an alternative
to the grouped $L_0/L_1$ representation of
[GenericInterstitialLoopGroupedNodalKernel.md]: the physics (loop radius, size-dependent bias,
detailed-balance emission from a capillary-law binding energy) is the same, and only the way the
large sizes are stored differs.

!alert warning title=One family only
This kernel has no opposite-species absorption and no recombination. It models a single loop
family with a monomer source and a linear sink.

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

!syntax parameters /NodalKernels/ContinuumLoopNodalKernel

!syntax inputs /NodalKernels/ContinuumLoopNodalKernel

!syntax children /NodalKernels/ContinuumLoopNodalKernel
