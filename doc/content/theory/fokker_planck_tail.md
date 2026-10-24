# Fokker-Planck Continuum Tail

The cluster size distributions of [DefectDynamics](syntax/DefectDynamics/index.md) and of
[SoluteClusterDynamics](syntax/SoluteClusterDynamics/index.md) can store the large cluster sizes in two ways. In the
default, the sizes above `explicit_max` are grouped in bins with a linear reconstruction ($L_0$, $L_1$). With
`tail = fokker_planck`, sizes up to `switch_size` are individual unknowns and the larger sizes are finite-volume
cells of a Fokker-Planck continuum, as described by [!cite](pokor2004). This page records what the option does,
what was measured with it, and what is not established. All numbers were measured in October 2026 on one laptop (17 GB of memory)
and are not general.

## The scheme

Sizes $n = 1, \ldots, n_s$ ($n_s$ is `switch_size`) are unknowns. Above $n_s$ the concentration per unit size is
averaged over cells that start at $n_s + 1/2$ with width `continuum_initial_width`, grow by `continuum_width_ratio`
(default 1.04) and end at $N + 1/2$, where $N$ is the largest size. The flux between two sizes is
$J_n = \beta_n C_1 C_n - \alpha_{n+1} C_{n+1}$. Between two cells it is an upwinded advective part, with drift
$\beta C_1 - \alpha$, and a diffusive part with the coefficient $(\beta C_1 + \alpha)/2$. The face values are reconstructed with a
minmod-limited slope (`continuum_limiter`, second order) or are the cell averages (first order). The flux out of the last cell is
zero, so $N$ must be large enough that the population does not reach it. When two species are present, each shrinks by absorbing a
monomer of the other, with its own upwinding, and the Jacobian, including the off-diagonal block, is written out by hand and
includes the dependence of the limited slopes on the neighboring cells.

The coefficients are those of the grouped kernels. Implemented: loops and voids with the capillary binding energy
([ContinuumLoopNodalKernel.md]); solute clusters with the `simple` and `interfacial_energy` rate models and the vacancy
coupling ([ContinuumClusterDynamicsNodalKernel.md]). Outputs: [ClusterContinuumMoment.md] and
[ClusterContinuumDistribution.md].

## What was verified

- The rate agrees with the SUNDIALS Pokor benchmark's right-hand side on the same state to $10^{-14}$ (monomers) and $3 \times 10^{-14}$ (cells);
  the Jacobian agrees with a finite difference to $2 \times 10^{-8}$ (PETSc test: $5 \times 10^{-9}$).
- The `[DefectDynamics]` and `[SoluteClusterDynamics]` options reproduce the same problems built from the kernels to $10^{-12}$ or better.
- Existing behavior is unchanged: the outputs of 51 existing inputs are byte-identical, and the existing tests pass.
- Against the SUNDIALS Pokor benchmark with identical physics (no grain-boundary sink, the 2/1 dimer bookkeeping), Fokker-Planck
  with `switch_size = 100` agrees to within 0.02 percent in the monomers, densities and mean radii at 0.5 to 30 dpa (30 dpa, largest size $6 \times 10^6$).

## Accuracy: the time step decides, not the scheme

For six void-and-loop cases (Yang, MD and Pokor parameters at 573 and 673 K, to 10 dpa, largest size 998,402) the grouped and
the Fokker-Planck runs were compared with the SUNDIALS geometric-grid reference. The physics parameters and equations of that
reference were checked to be those of BlackBear, but it is a grouped solution itself, so it is not ground truth.

- With the time step of the production inputs (`IterationAdaptiveDT` with a cap proportional to the time), **both** schemes are wrong
  in the resolvable densities during their steep rise, by factors of 2 to $10^5$ (rows where the reference density exceeds $10^{12}$ m$^{-3}$).
  An earlier result that Fokker-Planck was clearly more accurate came from this cap and is withdrawn.
- With the step also limited by the relative change of the resolvable densities (0.1 per step; recipe in the
  [DefectDynamics](syntax/DefectDynamics/index.md) page), the two are equally accurate. The median over the six cases of the worst deviation is
  12.1 percent (grouping) and 12.2 percent (Fokker-Planck, `switch_size = 160`) in the void density, 3.7 and 3.7 in the loop density,
  0.0 and 0.1 in the void radius, and 0.1 and 1.2 in the loop radius. The remaining void-density difference of about 11 percent is shared by both
  and lies in the steep onset; it was not traced to the reference or to the time discretization.
- For the Cu precipitation case (100,000 sizes, 563 K) the same holds: without a step cap both are off by a factor of 3 in the
  density at 0.01 s; with the step limited in proportion to the time, grouped and continuum agree with the SUNDIALS reference to
  within 0.4 and 1.3 percent in the density.
- `switch_size` must be at least the size at which the resolvable population starts (for a radius of 0.75 nm: 152 for voids, 32 for loops),
  otherwise the cut falls inside a cell and the resolvable void density is off by about 34 percent (0.3 percent at 160). The action warns.
  The cell ratio 1.02 against 1.04 changes the number densities by under 1 percent and doubles the memory.

## Cost: this is the advantage

When array variables are coupled, MOOSE allocates storage for every pair of components of the whole nonlinear system, so the memory
grows with the square of the number of unknowns. The continuum needs far fewer unknowns for the same accuracy.

| Problem | Grouping | Fokker-Planck |
|---|---|---|
| Six void+loop cases, to 10 dpa, adaptive step | 4268 unknowns, 523 to 729 s, 29.9 GB | 860 unknowns (switch 160), 10.6 to 15 s, 1.46 GB |
| The same, fixed step cap | 61 to 79 s | 6.4 to 7.0 s (switch 100) |
| Solute coupled to the defect arrays, production size | 4606 unknowns, 354.6 s, 35.3 GB | 1271 unknowns, 8.2 s, 2.98 GB |
| Cu, 100,000 sizes, solute alone | 338 unknowns, 0.04 GB, 4.3 s | 311 to 584 unknowns, 0.03 GB, about 2 s |

The memory is the **footprint** (dirty plus compressed pages), not the resident size: the grouped runs show a maximum resident size of 2.6 to
4.6 GB but a footprint of 29.9 GB, because macOS compresses pages that were allocated and not touched. The footprint matches the
estimate $10 \times 2 \times (\text{unknowns})^2 \times 80$ B of the pair storage (29.2 GB for 4268 unknowns). Solute clustering alone gains no memory
(the grouped array is already small); the gain comes when it shares a problem with large point-defect arrays.

## Not established

- Convergence in the cell ratio and in the switch size beyond the values tried (switch 50 to 400, ratio 1.02 to 1.08, mostly at one case).
- Binding-energy tables, helium and multiple solutes with the continuum (the tables are rejected by the actions).
- Other temperatures and parameter sets than those above; other machines; threading and MPI.
- An adaptive time-step rule built into the action; the recipe must be written in the input.
- The origin of the 11 percent void-density difference that both schemes share.

## See Also

- [DefectDynamics](syntax/DefectDynamics/index.md)
- [SoluteClusterDynamics](syntax/SoluteClusterDynamics/index.md)
- [ContinuumLoopNodalKernel.md]
- [ContinuumClusterDynamicsNodalKernel.md]
