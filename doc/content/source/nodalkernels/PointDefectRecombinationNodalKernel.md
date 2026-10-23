# PointDefectRecombinationNodalKernel

!syntax description /NodalKernels/PointDefectRecombinationNodalKernel

## Description

See [theory/defect_dynamics.md] for a full derivation connecting this
kernel's equation to the rest of the defect-dynamics model.

`PointDefectRecombinationNodalKernel` adds point-defect recombination loss to
one defect family's monomer component (array component 0), coupled to the
other family's own array variable through
[!param](/NodalKernels/PointDefectRecombinationNodalKernel/coupled_var). It is
the object that adds recombination between the two families of the
defect-dynamics model. The void and loop kernels
([GenericVacancyVoidGroupedNodalKernel.md],
[GenericInterstitialLoopGroupedNodalKernel.md]) couple the families only
through the opposite-species absorption described on their pages, and each
otherwise evolves its own grouped array variable.

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
The per-component names can be used instead, for example
`off_diag_row = 'vacancy_clusters_0 interstitial_clusters_0'` with
`off_diag_column = 'interstitial_clusters_0 vacancy_clusters_0'` and `full = false`.
`full = true` is the simpler choice for small problems.

!alert warning title=Memory grows with the square of the cluster-size count
Declaring the coupling between two grouped array variables, by either route, makes MOOSE
allocate storage for the coupling between every pair of their array components. Measured
with the opposite-species absorption active, 20,000 sizes per species took about 3.7 GB and
14 s for two steps with `full = true`, against about 0.12 GB and 1.7 s with `full = false` and
no cross-species block. Without the cross-species block, Newton's method uses an approximate
Jacobian for the coupling between the families: it still converges to the same residual, but
took more nonlinear iterations in the cases tried. For production-size problems, choose
between the full Jacobian and the memory.

## Example Input Syntax

!listing test/tests/defect_dynamics/point_defect_recombination_small.i

!syntax parameters /NodalKernels/PointDefectRecombinationNodalKernel

!syntax inputs /NodalKernels/PointDefectRecombinationNodalKernel

!syntax children /NodalKernels/PointDefectRecombinationNodalKernel

## See Also

- [theory/defect_dynamics.md]
- [GenericVacancyVoidGroupedNodalKernel.md]
- [GenericInterstitialLoopGroupedNodalKernel.md]
