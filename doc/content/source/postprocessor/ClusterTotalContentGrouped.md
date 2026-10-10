# ClusterTotalContentGrouped

!syntax description /Postprocessors/ClusterTotalContentGrouped

## Description

`ClusterTotalContentGrouped` computes the total point-defect content
represented by an explicit and grouped cluster distribution: the number of
point defects (monomers, vacancies, or interstitials) tied up in clusters of
every tracked size, rather than the number of clusters themselves
([ClusterTotalDensityGrouped.md] computes the latter):

!equation
S = \frac{1}{V_{at}}\sum_{n=n_\mathrm{min}}^N n\,C_n.

For a complete bin without reconstruction clipping, the sum over the bin
reduces to $\Delta(L_0\bar n + L_1\sigma_n^2)$: the cross term in
$\sum_n(n-\bar n)n$ vanishes by the definition of $\bar n$, leaving
$\sum_n(n-\bar n)n = \Delta\sigma_n^2$. Other bins are summed over their
reconstructed physical concentrations. When `atomic_volume` is provided in
m$^3$, the result is reported in `#/m^3`.

This postprocessor provides the moment ($\sum n C_n$) needed to check
[!citep](golubov2001)'s grouping-validity criterion (c), that a grouping
method must conserve the total number of point defects accumulated in
clusters, $dS/dt = J(1,t) + \sum_x J(x,t)$; see [theory/defect_dynamics.md]
for the flux definition $J_n$. [ClusterTotalDensityGrouped.md] provides the
companion moment ($\sum C_n$) for criterion (b), conservation of the total
number of clusters.

The grouping layout and nonnegative-reconstruction parameters must match those
used by [GenericClusterDynamicsGroupedNodalKernel.md].

## Example Input Syntax

!listing test/tests/cluster_dynamics/cluster_dynamics_grouped_conservation.i block=Postprocessors/total_cluster_content

!syntax parameters /Postprocessors/ClusterTotalContentGrouped

!syntax inputs /Postprocessors/ClusterTotalContentGrouped

!syntax children /Postprocessors/ClusterTotalContentGrouped

## See Also

- [ClusterTotalDensityGrouped.md]
- [theory/defect_dynamics.md]
