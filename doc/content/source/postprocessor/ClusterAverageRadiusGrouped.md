# ClusterAverageRadiusGrouped

!syntax description /Postprocessors/ClusterAverageRadiusGrouped

## Description

`ClusterAverageRadiusGrouped` computes the number-weighted mean radius from an
explicit and grouped cluster distribution. Only physical cluster sizes
$n \geq n_\mathrm{min}$ are included.

!equation
\bar r =
\frac{\displaystyle\sum_{n=n_\mathrm{min}}^N r_1 n^{p} C_n}
{\displaystyle\sum_{n=n_\mathrm{min}}^N C_n}.

The exponent $p$ is
[!param](/Postprocessors/ClusterAverageRadiusGrouped/radius_exponent), $1/3$
by default for three-dimensional clusters such as the spherical clusters of
[GenericClusterDynamicsGroupedNodalKernel.md] and the voids of
[GenericVacancyVoidGroupedNodalKernel.md]. Set it to $1/2$ for
two-dimensional clusters such as the interstitial loops of
[GenericInterstitialLoopGroupedNodalKernel.md].

For complete grouped bins, the denominator is evaluated from the stored zeroth
coefficient. Radius sums needed by the numerator are precomputed for each bin.
Partially included bins and bins requiring reconstruction limiting are
evaluated over their physical cluster sizes.

The grouping layout and nonnegative-reconstruction parameters must match those
used by [GenericClusterDynamicsGroupedNodalKernel.md].

## Example Input Syntax

!listing test/tests/cluster_dynamics/cluster_dynamics_grouped_small.i block=Postprocessors/avg_cluster_radius

!syntax parameters /Postprocessors/ClusterAverageRadiusGrouped

!syntax inputs /Postprocessors/ClusterAverageRadiusGrouped

!syntax children /Postprocessors/ClusterAverageRadiusGrouped
