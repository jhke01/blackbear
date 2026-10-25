# DefectRecoveryAction

!syntax description /DefectDynamics/Recovery/DefectRecoveryAction

## Description

The `DefectRecoveryAction` is a MOOSE action that describes the shared dislocation network
parameters (initial density, recovery rate constant, Burgers vector) used by every species'
[GenericDislocationRecoveryNodalKernel.md] in the enclosing `DefectDynamics` block, when this
optional sub-block is given. This is invoked using the `Recovery` sub-block within
[DefectDynamics](syntax/DefectDynamics/index.md). See that page for a description, example
usage, and parameters.

## Models

[!param](/DefectDynamics/Recovery/DefectRecoveryAction/model) selects how the density evolves.

- `pokor` (default): $\rho$ recovers by climb, $d\rho/dt = -K b^2 \rho^{3/2}$ ([!cite](pokor2004)), in
  closed form; needs [!param](/DefectDynamics/Recovery/DefectRecoveryAction/K).
- `cellular_network`: the cellular network dislocation density of [!cite](rest2006). Interstitial loops
  join the network at a rate set by the interstitial monomer concentration $c_i$, and the network loses
  line length by self-annihilation and at grain boundaries ([CellularNetworkDislocationDensity.md]).
  It needs [!param](/DefectDynamics/Recovery/DefectRecoveryAction/bias) (Rest's $Z_{iv}$) and uses
  [!param](/DefectDynamics/Recovery/DefectRecoveryAction/ca_crho),
  [!param](/DefectDynamics/Recovery/DefectRecoveryAction/poisson_ratio) and
  [!param](/DefectDynamics/Recovery/DefectRecoveryAction/grain_size); $K$ is not allowed. The loop
  species whose monomer drives the network is the only species with `geometry = loop`, or is named by
  [!param](/DefectDynamics/Recovery/DefectRecoveryAction/interstitial_species).

With [!param](/DefectDynamics/Recovery/DefectRecoveryAction/loop_transfer) = true (and `tail = fokker_planck`) the
network is not fed by the steady-state loop density of Rest's Eq. 18. The interstitial loops of the cluster
distribution themselves join the network at the rate $k = v_l / L$ ([LoopNetworkTransferNodalKernel.md]), and the
network density gains the line length they carry. Loops of at least
[!param](/DefectDynamics/Recovery/DefectRecoveryAction/transfer_minimum_size) join. This removes the double counting
of line length between the explicit loops and the network. The rate $k$ for every loop above a minimum size is a modelling choice
of this implementation, not of [!cite](rest2006).

!alert warning title=loop_transfer is experimental
In the SUNDIALS reference of this model (SA304 parameters, $Z_{iv} = 0.1$, $C_A C_\rho = 2$, 10 dpa) the network density
reaches the same value of about $1.4 \times 10^{15}$ m$^{-2}$ from initial densities of $10^{10}$ and $10^{15}$ m$^{-2}$ at 603 K, but at 673 K it does not:
it ends at $5.5 \times 10^{9}$ m$^{-2}$ from $10^{10}$ and at $6 \times 10^{13}$ m$^{-2}$ from $10^{15}$. At 673 K with $Z_{iv} = 0.1$ there is a critical starting density between
$3 \times 10^{11}$ and $6 \times 10^{11}$ m$^{-2}$: lower starts fall, higher starts reach a plateau of $6 \times 10^{13}$ m$^{-2}$. At 573 and 623 K the two starts tried
($10^{10}$ and $10^{15}$) end within a factor of 1.5. Other values of $Z_{iv}$ and other temperatures have not been mapped. The results depend on the transfer rate, which is an assumption of this implementation.
Do not treat the density it gives as a prediction.

!alert warning title=The bias is not known for steel
Rest gives $Z_{iv} = 5 \times 10^{-4}$ for UO$_2$. No value for steel is given in that paper. In the tests and in
the comparison with the SUNDIALS reference, 0.1 and 0.01 are used as round numbers, and the result depends strongly
on this choice (a steady-state density that changes by a factor of 10 or more).
