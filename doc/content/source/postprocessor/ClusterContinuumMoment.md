# ClusterContinuumMoment

!syntax description /Postprocessors/ClusterContinuumMoment

## Description

`ClusterContinuumMoment` computes a moment of the size distribution of a cluster family that is
discrete up to [!param](/Postprocessors/ClusterContinuumMoment/switch_size) and a Fokker-Planck
continuum above it, the representation of [ContinuumLoopNodalKernel.md]. The
[!param](/Postprocessors/ClusterContinuumMoment/quantity) parameter selects

- `density`, the number density $\sum_n C_n / V_{at}$ of clusters,
- `content`, the density of point defects they hold, $\sum_n n C_n / V_{at}$, or
- `average_radius`, the number-weighted mean radius $\sum_n r_n C_n / \sum_n C_n$ with
  $r_n = r_1 n^p$ ([!param](/Postprocessors/ClusterContinuumMoment/r1) and
  [!param](/Postprocessors/ClusterContinuumMoment/radius_exponent)).

Only sizes of at least [!param](/Postprocessors/ClusterContinuumMoment/n_minimum) are counted.
Choosing the size that corresponds to the smallest radius a microscope can resolve gives the
observable population, as in the SUNDIALS Pokor benchmark. A continuum cell holds a concentration
per unit size, so its contribution is the integral over the part of the cell above the lower edge
$n_\mathrm{min} - 1/2$ of the first counted size; a cell that straddles that edge is cut there.
Concentrations enter as $\max(C, 0)$, as in the kernel.

The layout parameters
([!param](/Postprocessors/ClusterContinuumMoment/num_cluster_sizes),
[!param](/Postprocessors/ClusterContinuumMoment/switch_size),
[!param](/Postprocessors/ClusterContinuumMoment/continuum_initial_width) and
[!param](/Postprocessors/ClusterContinuumMoment/continuum_width_ratio)) must be those of the
cluster variable.

## Example Input Syntax

In this example a variable holds a constant concentration of 2 per unit size on every discrete size
and every cell, so the moments are known exactly. The first three postprocessors count every size,
and the last three count only sizes of 11 and above.

!listing test/tests/continuum_tail/continuum_moments.i block=Postprocessors

!syntax parameters /Postprocessors/ClusterContinuumMoment

!syntax inputs /Postprocessors/ClusterContinuumMoment

!syntax children /Postprocessors/ClusterContinuumMoment
