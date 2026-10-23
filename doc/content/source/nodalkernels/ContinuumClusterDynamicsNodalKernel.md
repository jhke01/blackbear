# ContinuumClusterDynamicsNodalKernel

!syntax description /NodalKernels/ContinuumClusterDynamicsNodalKernel

## Description

`ContinuumClusterDynamicsNodalKernel` evolves the size distribution of a solute-clustering family (Cu precipitation, for
example) with a representation that is discrete for small sizes and a Fokker-Planck continuum for large ones, as
described for large clusters by [!cite](pokor2004). It is an alternative to the grouped $L_0/L_1$ representation of
[GenericClusterDynamicsGroupedNodalKernel.md]: the rate coefficients are the same, and only the way the large sizes are
stored differs. The [SoluteClusterDynamics](syntax/SoluteClusterDynamics/index.md) block builds it with
[!param](/SoluteClusterDynamics/SoluteClusterDynamicsAction/tail) `= fokker_planck`.

## Representation

The array variable has one component for the monomer, one for each size $2, \ldots, n_s$
([!param](/NodalKernels/ContinuumClusterDynamicsNodalKernel/switch_size)) and one for each finite-volume cell above it
(the cells start at $n_s + 1/2$ with the width
[!param](/NodalKernels/ContinuumClusterDynamicsNodalKernel/continuum_initial_width), grow by
[!param](/NodalKernels/ContinuumClusterDynamicsNodalKernel/continuum_width_ratio) and end at
[!param](/NodalKernels/ContinuumClusterDynamicsNodalKernel/num_cluster_sizes) $+ 1/2$). The flux between two sizes is
$J_n = \beta_n C_1 C_n - \alpha_{n+1} C_{n+1}$; between two cells it is an advective part, upwinded, and a diffusive part,
with the coefficients interpolated linearly in the size, and
[!param](/NodalKernels/ContinuumClusterDynamicsNodalKernel/use_limiter) selects a second-order (minmod-limited) or a
first-order reconstruction of the face values. The flux out of the last cell is zero, so
[!param](/NodalKernels/ContinuumClusterDynamicsNodalKernel/num_cluster_sizes) has to be large enough that the population
does not reach it. The Jacobian is written out by hand and includes the dependence of the limited slopes on the
neighboring cells.

## Rate Coefficients

With [!param](/NodalKernels/ContinuumClusterDynamicsNodalKernel/rate_model) `= interfacial_energy`, clusters are spheres
of radius $r_n = (3 n V_{at} / 4 \pi)^{1/3}$, the absorption coefficient is
$\beta_n = 4 \pi (r_1 + r_n) D / V_{at}$ and the emission coefficient follows from detailed balance,
$\alpha_n = \beta_{n-1} \exp(-E_b(n) / k_B T)$, with the binding energy $E_b(n)$ from the interfacial energy
([!param](/NodalKernels/ContinuumClusterDynamicsNodalKernel/sigma),
[!param](/NodalKernels/ContinuumClusterDynamicsNodalKernel/Omega_kB_K),
[!param](/NodalKernels/ContinuumClusterDynamicsNodalKernel/DeltaS_kB)) or from the capillary law
(`binding_energy_model = capillary`). With `rate_model = simple`, $\beta_n = \beta_0 n^{1/3}$ and
$\alpha_n = \alpha_0 n^{1/3}$. The binding-energy tables of the grouped kernel are not supported.

## Vacancy Coupling

With [!param](/NodalKernels/ContinuumClusterDynamicsNodalKernel/vacancy_coupled_var), the monomer diffusivity is
multiplied by the vacancy supersaturation, the monomer concentration of the coupled array over its thermal value
$4 \exp(-E_f / k_B T)$ ([!param](/NodalKernels/ContinuumClusterDynamicsNodalKernel/vacancy_formation_energy_eV)). The
rates are proportional to the diffusivity, so the off-diagonal Jacobian column of the vacancy monomer is the rate without
the sources divided by that monomer. The preconditioner needs the pair declared (`off_diag_row`, `off_diag_column`).

## Example Input Syntax

In this example the solute array is discrete up to size 20 and a continuum up to size 2000, and its monomer diffusivity
is coupled to the vacancy monomer of a second array.

!listing test/tests/continuum_tail/continuum_solute_vacancy_jacobian.i block=NodalKernels/dyn

!syntax parameters /NodalKernels/ContinuumClusterDynamicsNodalKernel

!syntax inputs /NodalKernels/ContinuumClusterDynamicsNodalKernel

!syntax children /NodalKernels/ContinuumClusterDynamicsNodalKernel
