# ClusterSizeDistributionGrouped

!syntax description /VectorPostprocessors/ClusterSizeDistributionGrouped

## Description

`ClusterSizeDistributionGrouped` reports the full per-size cluster concentration distribution of
an explicit and grouped cluster array variable: one row (`n`, `concentration`, `radius`) for
every physical cluster size $n$ from 1 to
[!param](/VectorPostprocessors/ClusterSizeDistributionGrouped/num_cluster_sizes).

Explicit sizes ($n \leq$
[!param](/VectorPostprocessors/ClusterSizeDistributionGrouped/explicit_max)) are read directly
from their own array component. Grouped (large-$n$) sizes are reconstructed from their bin's
stored $L_0$/$L_1$ moments using the same linear fit
$C_n = a_0 L_0 + a_1 L_1$ the grouped kernels themselves use internally
(`ClusterGroupingLayout::concentration()`), optionally clipped so bin endpoint concentrations
stay nonnegative. Because every $n$ inside one bin shares the same $L_0$/$L_1$, adjacent
reconstructed rows within a bin lie on that bin's single line rather than being independent
data -- this object makes the full array readable at every physical size, not finer-grained than
what the grouping actually stores.

`radius` is $r_1 n^p$, with $r_1$ =
[!param](/VectorPostprocessors/ClusterSizeDistributionGrouped/r1) and $p$ =
[!param](/VectorPostprocessors/ClusterSizeDistributionGrouped/radius_exponent) ($1/3$ for a
three-dimensional void, $1/2$ for a two-dimensional loop).

The grouping layout and nonnegative-reconstruction parameters must match those used by the
coupled array variable's own kernel (e.g. [GenericVacancyVoidGroupedNodalKernel.md] or
[GenericInterstitialLoopGroupedNodalKernel.md]).

Since this model is run on a single-node (0D) mesh in every current production and test input,
one node's `execute()` call already covers the entire distribution; `finalize()` and
`threadJoin()` are no-ops.

In `[DefectDynamics]` syntax, a species adds this object by setting
[!param](/DefectDynamics/DefectSpeciesAction/output_size_distribution) `= true` -- see
[syntax/DefectDynamics/index.md#full-cluster-size-distribution] for the full input idiom and the
dedicated sync-times-only output it is written through.

## Example Input Syntax

!listing test/tests/defect_dynamics/defect_dynamics_size_distribution_action.i block=DefectDynamics

!syntax parameters /VectorPostprocessors/ClusterSizeDistributionGrouped

!syntax inputs /VectorPostprocessors/ClusterSizeDistributionGrouped

!syntax children /VectorPostprocessors/ClusterSizeDistributionGrouped

## See Also

- [ClusterTotalDensityGrouped.md]
- [ClusterAverageRadiusGrouped.md]
- [syntax/DefectDynamics/index.md]
- [theory/defect_dynamics.md]
