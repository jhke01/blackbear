# ClusterContinuumDistribution

!syntax description /VectorPostprocessors/ClusterContinuumDistribution

## Description

`ClusterContinuumDistribution` reports the size distribution of a cluster family that is discrete
up to [!param](/VectorPostprocessors/ClusterContinuumDistribution/switch_size) and a Fokker-Planck
continuum above it, the representation of [ContinuumLoopNodalKernel.md]. It has one row for every
discrete size and one for every continuum cell, with the columns

- `n`, the size, or the center of the size interval for a cell,
- `width`, the width of that interval (1 for a discrete size),
- `concentration`, the concentration per unit size, and
- `radius`, $r_1 n^p$ at `n` ([!param](/VectorPostprocessors/ClusterContinuumDistribution/r1) and
  [!param](/VectorPostprocessors/ClusterContinuumDistribution/radius_exponent)).

Multiplying `concentration` by `width` gives the number of clusters in the interval. Unlike
[ClusterSizeDistributionGrouped.md], which reports every physical size, the number of rows is the
number of unknowns, so the output stays small for a ceiling of millions of atoms.

The layout parameters must be those of the cluster variable.

## Example Input Syntax

In this example the variable holds a constant concentration of 2 per unit size, with sizes up to
20 discrete and a continuum up to 2000.

!listing test/tests/continuum_tail/continuum_distribution.i block=VectorPostprocessors

!syntax parameters /VectorPostprocessors/ClusterContinuumDistribution

!syntax inputs /VectorPostprocessors/ClusterContinuumDistribution

!syntax children /VectorPostprocessors/ClusterContinuumDistribution
