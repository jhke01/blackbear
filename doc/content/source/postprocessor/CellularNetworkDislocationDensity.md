# CellularNetworkDislocationDensity

!syntax description /Postprocessors/CellularNetworkDislocationDensity

## Description

`CellularNetworkDislocationDensity` advances the dislocation density $\rho_N$ of a cellular dislocation network
by the model of [!cite](rest2006), Eqs. 14 to 18. The network gains line length from interstitial loops
and loses it by self-annihilation of adjacent edge dislocations and at grain boundaries:

\begin{equation}
\frac{d\rho_N}{dt} = \pi v_l n_l - \frac{4 v_l}{d_g}\rho_N - v_l \frac{(f(\nu)/\pi)^{1/2}}{C_A C_\rho}\rho_N^{3/2},
\qquad v_l = \frac{2 Z_{iv} D_i c_i}{b}, \quad f(\nu) = \frac{1 - \nu/2}{1 - \nu},
\end{equation}

where $n_l$ is the interstitial loop density, $D_i c_i$ the interstitial diffusivity and site fraction,
$Z_{iv}$ the relative bias ([!param](/Postprocessors/CellularNetworkDislocationDensity/bias)), $b$ the Burgers vector
([!param](/Postprocessors/CellularNetworkDislocationDensity/burgers_vector)), $d_g$ the grain size
([!param](/Postprocessors/CellularNetworkDislocationDensity/grain_size); the term is dropped for 0),
$\nu$ the Poisson ratio ([!param](/Postprocessors/CellularNetworkDislocationDensity/poisson_ratio)) and $C_A C_\rho$
the cell geometry factor ([!param](/Postprocessors/CellularNetworkDislocationDensity/ca_crho)).

The loop density is taken at its steady state, Rest's Eq. 18,

\begin{equation}
n_l = \left[\frac{\sqrt{2}\, b\, c_i\, \rho_N}{4\pi Z_{iv}\Omega^{5/3}}\right]^{1/2},
\end{equation}

with $\Omega$ the atomic volume ([!param](/Postprocessors/CellularNetworkDislocationDensity/atomic_volume)). With
$u = \sqrt{\rho_N}$ the equation is the Riccati equation $du/dt = c - b u - a u^2$, which is solved exactly over a step for a constant $c_i$.
The postprocessor uses the mean of $c_i$ ([!param](/Postprocessors/CellularNetworkDislocationDensity/interstitial_concentration))
at the two ends of the step. It cannot become negative. The accuracy is the accuracy of holding $c_i$ constant, so the
time step has to follow the change of $c_i$ and of $\rho_N$ (`timestep_relative_change` of the
[DefectDynamics](syntax/DefectDynamics/index.md) block includes the density).

With [!param](/Postprocessors/CellularNetworkDislocationDensity/line_length_postprocessor) the source is not the steady-state loop density but the line
length $L$ of the loops that join the network, and

\begin{equation}
\frac{d\rho_N}{dt} = k\,L - k\,\rho_N - \frac{4 v_l}{d_g}\rho_N, \qquad k = v_l \frac{\sqrt{f/\pi}}{C_A C_\rho}\sqrt{\rho_N},
\end{equation}

where $k\rho_N$ is Rest's self-annihilation term, so this is the same Riccati equation with the source $(C v_l / 2) L$ in the variable $\sqrt{\rho_N}$. [LoopNetworkTransferNodalKernel.md] removes the loops
from the cluster distribution at the same rate $k$.

The steady-state loop density is a good approximation late in the irradiation but not early: in the
SUNDIALS reference with the loop density as its own equation (Rest's Eq. 15), the network density at 0.1 dpa is about 25 percent lower than with
the steady state, and the two agree to 4 percent at 10 dpa.

!alert warning title=Parameters
$Z_{iv}$ and $C_A C_\rho$ are material and geometry parameters that are not known for steel. The value of
$Z_{iv}$ for UO$_2$ in the paper is $5 \times 10^{-4}$.

## Example Input Syntax

In this example the postprocessor is built by the `[Recovery]` sub-block of a `[DefectDynamics]` block with `model = cellular_network`,
and the dislocation sink of both species reads it.

!listing test/tests/defect_dynamics/defect_dynamics_cellular_network_action.i block=DefectDynamics

!syntax parameters /Postprocessors/CellularNetworkDislocationDensity

!syntax inputs /Postprocessors/CellularNetworkDislocationDensity

!syntax children /Postprocessors/CellularNetworkDislocationDensity
