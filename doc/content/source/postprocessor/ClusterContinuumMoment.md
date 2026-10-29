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
  [!param](/Postprocessors/ClusterContinuumMoment/radius_exponent)),
- `line_length`, the dislocation line length per volume of a loop family, $\sum_n 2\pi r_n C_n / V_{at}$
  in m$^{-2}$ ([!cite](rest2006), Eq. 9, $\rho = \pi n_l d_l$ with $d_l$ the loop diameter), or
- `lattice_strain`, the relative change of the lattice parameter caused by the loops,
  $\frac{1}{3}\sum_n n C_n$ ([!cite](rest2006), Eq. 4, $\pi b n_l d_l^2 / 12$). For a loop of $n$ atoms,
  $r_n^2 = n V_{at} / (\pi b)$, so $\pi b\, r_n^2 = n V_{at}$ and Rest's expression is one third of the
  volume fraction of the atoms in loops.

`line_length` and `lattice_strain` are defined for loops only
([!param](/Postprocessors/ClusterContinuumMoment/radius_exponent) = 0.5), and the monomer is not a
loop, so they are meant to be used with [!param](/Postprocessors/ClusterContinuumMoment/n_minimum) = 2.
The lattice strain is the relative change of the lattice parameter of the crystal that contains the loops, from the loop population only; it ignores other
sources of strain, and it assumes that each loop atom adds its atomic volume to the crystal.

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
