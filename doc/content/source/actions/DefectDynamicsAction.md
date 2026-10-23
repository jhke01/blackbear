# DefectDynamicsAction

!syntax description /DefectDynamics/DefectDynamicsAction

## Description

The `DefectDynamicsAction` is a MOOSE action that builds the grouped array variables, initial
conditions, self/opposite-absorption kernels, point-defect recombination kernels, and
monomer/TEM-visible postprocessors for one or two defect species, deriving every quantity the
underlying kernels need (sink coefficients, recombination rate, capillary coefficients, monomer
initial condition, TEM-visible cutoffs, and the opposite-species coupling parameters) from
physical inputs instead of requiring them to be typed and kept consistent by hand. This is
invoked using the [DefectDynamics](syntax/DefectDynamics/index.md) block. See that page for a
description, example usage, and parameters. With `tail = fokker_planck` the large cluster sizes are
a Fokker-Planck continuum ([ContinuumLoopNodalKernel.md]) instead of grouped bins.
