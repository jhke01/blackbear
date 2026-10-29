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

!alert warning title=Experimental
In the SUNDIALS reference of this model (SA304 parameters, $Z_{iv} = 0.1$, $C_A C_\rho = 2$, 10 dpa) the network density
reaches the same value of about $1.4 \times 10^{15}$ m$^{-2}$ from initial densities of $10^{10}$ and $10^{15}$ m$^{-2}$ at 603 K, but at 673 K it does not:
it ends at $5.5 \times 10^{9}$ m$^{-2}$ from $10^{10}$ and at $6 \times 10^{13}$ m$^{-2}$ from $10^{15}$. At 673 K with $Z_{iv} = 0.1$ there is a critical starting density between
$3 \times 10^{11}$ and $6 \times 10^{11}$ m$^{-2}$: lower starts fall, higher starts reach a plateau of $6 \times 10^{13}$ m$^{-2}$. At 573 and 623 K the two starts tried
($10^{10}$ and $10^{15}$) end within a factor of 1.6 (1.53 at 573 K, 1.07 at 623 K). Other values of $Z_{iv}$ and other temperatures have not been mapped. The results depend on the transfer rate, which is an assumption of this implementation.
Do not treat the density it gives as a prediction.

The kernel is added by `[Recovery]` with `loop_transfer = true` (see [DefectRecoveryAction.md]); it needs `tail = fokker_planck`.

## Example Input Syntax

In this example the loop species of a `[DefectDynamics]` block loses loops to the network. The kernel and the
line-length postprocessor are built by the `[Recovery]` sub-block.

!listing test/tests/defect_dynamics/defect_dynamics_loop_transfer_action.i block=DefectDynamics/Recovery

!syntax parameters /NodalKernels/LoopNetworkTransferNodalKernel

!syntax inputs /NodalKernels/LoopNetworkTransferNodalKernel

!syntax children /NodalKernels/LoopNetworkTransferNodalKernel
