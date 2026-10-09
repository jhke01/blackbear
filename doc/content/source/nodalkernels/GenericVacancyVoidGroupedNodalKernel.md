# ADVacancyVoidGroupedNodalKernel / VacancyVoidGroupedNodalKernel

!syntax description /NodalKernels/VacancyVoidGroupedNodalKernel

## Description

`VacancyVoidGroupedNodalKernel` solves the self-reaction (absorption and thermal
emission) terms of a vacancy-void population balance, using the same grouped
L0/L1 representation as [GenericClusterDynamicsGroupedNodalKernel.md]: cluster
sizes through
[!param](/NodalKernels/VacancyVoidGroupedNodalKernel/explicit_max) are stored
explicitly, and the remaining physical sizes are divided into bins carrying two
moments each. See that page for the grouped representation, moment evolution
equations, reconstruction limiter, and intra-variable Jacobian, which this
kernel reuses unchanged.

This kernel implements only the vacancy family's own growth/shrinkage by
vacancy-monomer absorption and emission. Point-defect recombination with an
interstitial variable and dislocation sink loss for the monomer (via the
[!param](/NodalKernels/VacancyVoidGroupedNodalKernel/sink) parameter) are the
only inter-family or sink physics present; interstitial-cluster capture by
voids is not yet implemented.

## Rate Coefficients

Voids are three-dimensional clusters, so the void radius follows

!equation
r_n = \left(\frac{3nV_{at}}{4\pi}\right)^{1/3},

the same functional form as the spherical-cluster radius in
[GenericClusterDynamicsNodalKernel.md]. The vacancy-monomer diffusivity follows
an Arrhenius law,

!equation
D_v = D_{0,v}\exp\left(-\frac{E_{m,v}}{k_BT}\right),

and the absorption rate coefficient includes a finite-sink-strength
correction that does not appear in the Cu-precipitation kernel:

!equation
\beta_n = \frac{4\pi r_n D_v}{\left(1+\dfrac{a_0}{r_n}\right)V_{at}},

where $a_0$ is
[!param](/NodalKernels/VacancyVoidGroupedNodalKernel/lattice_parameter). This
correction becomes significant once $r_n$ approaches $a_0$ and suppresses the
absorption rate relative to the uncorrected $4\pi r_n D_v/V_{at}$ form.

Thermal emission follows detailed balance from a capillary-law formation
energy,

!equation
E_f(n) = A_v n^{p_v} + C_v,
\qquad
E_b(n) = E_{form,v} - \left[E_f(n)-E_f(n-1)\right],
\qquad
\alpha_n = \beta_{n-1}\exp\left(-\frac{E_b(n)}{k_BT}\right),

where $p_v$ is
[!param](/NodalKernels/VacancyVoidGroupedNodalKernel/void_energy_exponent),
$A_v$ and $C_v$ are
[!param](/NodalKernels/VacancyVoidGroupedNodalKernel/void_energy_coefficient_eV)
and
[!param](/NodalKernels/VacancyVoidGroupedNodalKernel/void_energy_constant_eV),
and $E_{form,v}$ is
[!param](/NodalKernels/VacancyVoidGroupedNodalKernel/formation_energy_eV).

## Required and Recommended Solver Settings

This kernel requires `use_hash_table_matrix_assembly = true` in the
`[Problem]` block, for the same reason as
[GenericClusterDynamicsGroupedNodalKernel.md]. See that page for the full
discussion of recommended `[Preconditioning]` and `[Executioner]` settings.

## Example Input Syntax

!listing test/tests/defect_dynamics/vacancy_void_grouped_small.i

!syntax parameters /NodalKernels/VacancyVoidGroupedNodalKernel

!syntax inputs /NodalKernels/VacancyVoidGroupedNodalKernel

!syntax children /NodalKernels/VacancyVoidGroupedNodalKernel

## See Also

- [GenericClusterDynamicsGroupedNodalKernel.md]
- [GenericInterstitialLoopGroupedNodalKernel.md]
- [GroupedArrayFunctionIndexIC.md]
