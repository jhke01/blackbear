# ADDislocationRecoveryNodalKernel / DislocationRecoveryNodalKernel

!syntax description /NodalKernels/DislocationRecoveryNodalKernel

## Description

See [theory/defect_dynamics.md] for a full derivation connecting this kernel's equation to the
rest of the defect-dynamics model.

`DislocationRecoveryNodalKernel` adds dislocation-sink loss to one defect species' monomer
component (array component 0), using a dislocation density that evolves in time via
irradiation-induced climb instead of the fixed value the core void/loop kernels' own
[!param](/NodalKernels/VacancyVoidGroupedNodalKernel/sink) parameter assumes. It replaces that
constant sink term for species where this module is active -- the two are not meant to be used
together on the same species (see [DefectDynamics](syntax/DefectDynamics/index.md), which
enforces this at the action level).

## Equation

Pokor et al. (2004) attribute the vanishing of the initial dislocation network under irradiation
to dislocation climb driven by net absorption of point defects. The dislocation density $\rho$
then evolves as

!equation
\frac{d\rho}{dt} = -K b^2 \rho^{3/2},

with the closed-form solution

!equation
\rho(t) = \left(\rho_0^{-1/2} + \frac{1}{2}K b^2 t\right)^{-2},

where $\rho_0$ is [!param](/NodalKernels/DislocationRecoveryNodalKernel/rho0), $K$ is
[!param](/NodalKernels/DislocationRecoveryNodalKernel/K), and $b$ is
[!param](/NodalKernels/DislocationRecoveryNodalKernel/burgers_vector) -- the dislocation
network's own Burgers vector, a distinct quantity from any loop species' own geometric Burgers
vector (e.g. [!param](/NodalKernels/InterstitialLoopGroupedNodalKernel/burgers_vector)). Neither
$\rho$, $K$, nor this $b$ carries a species subscript in Pokor's own equations: one shared
dislocation network feeds both species' sink terms, each through that species' own diffusivity
and bias factor,

!equation
\text{loss} = D Z \rho(t) C_1,

matching the structure of [GenericVacancyVoidGroupedNodalKernel.md]/
[GenericInterstitialLoopGroupedNodalKernel.md]'s own constant-sink term (`sink =
D*Z*dislocation_density`) but with $\rho(t)$ in place of a fixed density. Because
$\rho(t)$ depends only on the current simulation time, this kernel needs no coupling to any
other variable and no tracked state beyond $\rho_0$, $K$, and $b$ -- unlike
[PointDefectRecombinationNodalKernel.md], it is AD-safe and provided in both AD and non-AD
variants.

!alert note
Pokor's own text states $K$ in units of m$^{-1}$ for the three reactor conditions it quotes (300,
970, 14400), but $d\rho/dt = -Kb^2\rho^{3/2}$ requires $K$ in m$^{-1}$s$^{-1}$ for the equation to
balance dimensionally. Treat [!param](/NodalKernels/DislocationRecoveryNodalKernel/K) as a fit
parameter in self-consistent SI units, not a literal copy of the paper's quoted numbers.

## Example Input Syntax

!listing test/tests/defect_dynamics/dislocation_recovery_small.i block=NodalKernels/vacancy_recovery

!syntax parameters /NodalKernels/DislocationRecoveryNodalKernel

!syntax inputs /NodalKernels/DislocationRecoveryNodalKernel

!syntax children /NodalKernels/DislocationRecoveryNodalKernel

## See Also

- [theory/defect_dynamics.md]
- [GenericVacancyVoidGroupedNodalKernel.md]
- [GenericInterstitialLoopGroupedNodalKernel.md]
- [PointDefectRecombinationNodalKernel.md]
