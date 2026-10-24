# LoopNetworkTransferNodalKernel

!syntax description /NodalKernels/LoopNetworkTransferNodalKernel

## Description

`LoopNetworkTransferNodalKernel` removes interstitial loops from a cluster size distribution at the rate at
which they join the dislocation network. It is the counterpart of the line-length source of
[CellularNetworkDislocationDensity.md] when that postprocessor has a
[!param](/Postprocessors/CellularNetworkDislocationDensity/line_length_postprocessor): the loops that leave
the distribution here are the line length that the network gains there, so the line length leaves the loops and enters the network. The network then loses line length only by
self-annihilation (the $\rho_N^{3/2}$ term of Rest's Eq. 14) and at grain boundaries.

A loop climbs with the velocity $v_l = 2 Z_{iv} D_i c_i / b$ ([!cite](rest2006), Eq. 16), where $c_i$ is the monomer
(array component 0) of the same variable, and reaches the network after the distance
$L = C_A C_\rho \sqrt{\pi / (\rho f(\nu))}$ (the cell size, Eq. 10) to the nearest wall. Every loop of at least
[!param](/NodalKernels/LoopNetworkTransferNodalKernel/minimum_size) therefore leaves its size at the rate

\begin{equation}
k = \frac{v_l}{L} = v_l \frac{\sqrt{f(\nu)/\pi}}{C_A C_\rho}\sqrt{\rho},
\end{equation}

the same for every size, so the residual is $k C_n$ on every component from the minimum size up (discrete sizes and
continuum cells). The density $\rho$ is read from [!param](/NodalKernels/LoopNetworkTransferNodalKernel/dislocation_density_postprocessor)
and is constant within a step; the dependence of $k$ on the monomer is part of the Jacobian.

!alert warning title=A modelling choice, not a result of the paper
Rest's model has one loop population with a diameter that is tied to the cell size (his Eq. 9 and 10) and does not
say how loops of the many sizes of a cluster-dynamics distribution join the network. The rate $k = v_l / L$ for every loop
above a minimum size is an assumption made here. Loops are moved to the network at a rate that does not depend on their own size or on
the distance of that particular loop to a wall.

The kernel is added by `[Recovery]` with `loop_transfer = true` (see [DefectRecoveryAction.md]); it needs `tail = fokker_planck`.

## Example Input Syntax

In this example the loop species of a `[DefectDynamics]` block loses loops to the network. The kernel and the
line-length postprocessor are built by the `[Recovery]` sub-block.

!listing test/tests/defect_dynamics/defect_dynamics_loop_transfer_action.i block=DefectDynamics/Recovery

!syntax parameters /NodalKernels/LoopNetworkTransferNodalKernel

!syntax inputs /NodalKernels/LoopNetworkTransferNodalKernel

!syntax children /NodalKernels/LoopNetworkTransferNodalKernel
