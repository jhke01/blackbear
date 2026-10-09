# PointDefectRecombinationNodalKernel

!syntax description /NodalKernels/PointDefectRecombinationNodalKernel

## Description

`PointDefectRecombinationNodalKernel` adds point-defect recombination loss to
one defect family's monomer component (array component 0), coupled to the
other family's own array variable through
[!param](/NodalKernels/PointDefectRecombinationNodalKernel/coupled_var). It is
the only object in the defect-dynamics model
([GenericVacancyVoidGroupedNodalKernel.md],
[GenericInterstitialLoopGroupedNodalKernel.md]) that couples the two families
to each other; each kernel otherwise evolves only its own grouped array
variable.

Recombination consumes one vacancy monomer and one interstitial monomer per
event, at rate

!equation
R = k_{rec} C_{v,1} C_{i,1},
\qquad
k_{rec} = \frac{4\pi r_{rec}(D_i+D_v)}{V_{at}},

matching the reference AMSS defect-dynamics model
(`DefectDynamicsAMSSInput`/`DefectDynamics_Sundials.cpp` in
`CD_data_for_comparison_mathew`). Two symmetric instances are required to
represent recombination: one applied to the vacancy variable (coupled to the
interstitial variable) and one applied to the interstitial variable (coupled
to the vacancy variable), both using the same
[!param](/NodalKernels/PointDefectRecombinationNodalKernel/recombination_rate).

## Non-AD Only

Unlike [GenericVacancyVoidGroupedNodalKernel.md] and
[GenericInterstitialLoopGroupedNodalKernel.md], this kernel does not provide
an AD variant. MOOSE's only array-variable coupling accessor,
`Coupleable::coupledArrayValue()`, returns a plain `ArrayVariableValue` with
no AD derivative tracking, and `ADArrayNodalKernel::setJacobian()`
unconditionally errors. An AD version of this kernel could therefore never
actually assemble its off-diagonal Jacobian block (the derivative with
respect to the *coupled* variable): it would silently leave that block zero,
degrading Newton's convergence rate without any compile- or run-time
indication that something was wrong. The residual, and therefore the
converged solution, is unaffected by this limitation; only automatic
differentiation with respect to the coupled variable is unavailable.

## Required Preconditioning

Because this kernel genuinely couples two different nonlinear variables, the
preconditioner's coupling matrix must include the inter-variable block.
`[Preconditioning]/[SMP]`'s `off_diag_row`/`off_diag_column` and
`coupled_groups` parameters resolve variable names through
`libMesh::System::variable_number()`, which does not recognize an array
variable's logical name (only its per-component names, e.g.
`vacancy_clusters_0`) -- a MOOSE-side limitation unrelated to this kernel.
Until that is addressed upstream, `full = true` is the practical way to
obtain a complete coupling matrix for two grouped array variables; see
[GenericClusterDynamicsGroupedNodalKernel.md] for why `full = false` matters
at production cluster-size counts for a single grouped variable, a
consideration this kernel's two-variable coupling does not change.

## Example Input Syntax

!listing test/tests/defect_dynamics/point_defect_recombination_small.i

!syntax parameters /NodalKernels/PointDefectRecombinationNodalKernel

!syntax inputs /NodalKernels/PointDefectRecombinationNodalKernel

!syntax children /NodalKernels/PointDefectRecombinationNodalKernel

## See Also

- [GenericVacancyVoidGroupedNodalKernel.md]
- [GenericInterstitialLoopGroupedNodalKernel.md]
