# ADClusterDynamicsNodalKernel / ClusterDynamicsNodalKernel

!syntax description /NodalKernels/ClusterDynamicsNodalKernel

## Description

The `ClusterDynamicsNodalKernel` implements the complete cluster dynamics rate equations for all cluster sizes (1 through N) in a +single array variable+. The rate coefficients can be supplied by either a simple power-law model or an interfacial-energy model based on cluster geometry and detailed balance for a single diffusing species.

!alert note title=Nodal Array Indexing
Array component index $i$ corresponds to cluster size $n = i+1$, i.e. monomer ($n=1$) array index is $i=0$; dimer ($n=2$) array index is $i=1$.  Largest cluster ($n=N$) array index is $i=N-1$.

### Rate Equations

+Component 0 (monomer, $n=1$):+

!equation
\frac{dC_1}{dt} = G_1 - k_s C_1 - 2\beta_1 C_1^2 - \sum_{n=2}^{N-1} \beta_n C_1 C_n + 2\alpha_2 C_2 + \sum_{n=3}^{N} \alpha_n C_n

where $G_1$ ([!param](/NodalKernels/ClusterDynamicsNodalKernel/generation)) is the monomer generation rate,  $k_s$ ([!param](/NodalKernels/ClusterDynamicsNodalKernel/sink)) is the  Linear sink coefficient

+Components+ $i \geq 1$ +(cluster size+ $n = i+1 \geq 2$+):+

!equation
\frac{dC_n}{dt} = \dot{C}_n^{\text{in}} - \beta_n C_1 C_n + \alpha_{n+1} C_{n+1} - \alpha_n C_n

where the growth-in term is:

!equation
\dot{C}_n^{\text{in}} = \beta_{n-1} C_1 C_{n-1}, \qquad n \geq 2

For $n = 2$, this reduces to $\beta_1 C_1^2$ because $C_{n-1} = C_1$.

!alert note title=Largest Cluster Truncation
For the largest tracked cluster size $n = N$, the system is closed at the upper bound of the
truncated cluster space. In that case, the forward absorption term to an untracked $N+1$ cluster
is set to zero and there is no emission-in term from $N+1$:

!equation
\frac{dC_N}{dt} = \beta_{N-1} C_1 C_{N-1} - \alpha_N C_N

This avoids an unphysical loss of mass from the tracked system through the top cluster bin and
preserves mass within the truncated cluster space.

The rate coefficients depend on the selected [!param](/NodalKernels/ClusterDynamicsNodalKernel/rate_model).

!equation
\beta_n = \beta_0 n^{1/3}, \qquad \alpha_n = \alpha_0 n^{1/3}

For `rate_model = simple`, $\beta_0$ ([!param](/NodalKernels/ClusterDynamicsNodalKernel/beta0)) and $\alpha_0$ ([!param](/NodalKernels/ClusterDynamicsNodalKernel/alpha0)) are user-supplied base coefficients.

For `rate_model = interfacial_energy`, the absorption coefficient is computed from the cluster geometry and monomer diffusivity:

!equation
V_{at} = \text{atomic volume}, \qquad r_n = \left(\frac{3 n V_{at}}{4\pi}\right)^{1/3}

!equation
\beta_n = \frac{4\pi (r_1 + r_n) D_m}{V_{at}}

where $V_{at}$ ([!param](/NodalKernels/ClusterDynamicsNodalKernel/atomic_volume)) is the atomic volume and $D_m$ is the monomer diffusivity. The diffusivity is selected by [!param](/NodalKernels/ClusterDynamicsNodalKernel/diffusivity_model), which is used only when `rate_model = interfacial_energy`:

!equation
D_m = \text{monomer diffusivity}

for `diffusivity_model = constant`, where $D_m$ is supplied directly by [!param](/NodalKernels/ClusterDynamicsNodalKernel/monomer_diffusivity).

!equation
D_m = D_0 \exp\left(-\frac{Q}{k_B T}\right)

for `diffusivity_model = arrhenius`, where $D_0$ ([!param](/NodalKernels/ClusterDynamicsNodalKernel/D0)) is the diffusion prefactor and $Q$ ([!param](/NodalKernels/ClusterDynamicsNodalKernel/Q)) is the activation energy.

In both diffusivity modes, the computed monomer diffusivity is multiplied by
[!param](/NodalKernels/ClusterDynamicsNodalKernel/radiation_enhanced_factor), which defaults to
1.0 and may be used to represent radiation-enhanced transport without manually rescaling the input
diffusivity or diffusion prefactor.

The emission coefficient is then derived by detailed balance using the cluster binding energy:

!equation
\alpha_n = \beta_{n-1}\exp\left(-\frac{E_n^b}{k_B T}\right), \qquad n \geq 2

!equation
E_n^b = \Omega - T\Delta S - (36\pi)^{1/3}V_{at}^{2/3}\sigma\left[n^{2/3} - (n-1)^{2/3}\right]

where $T$ ([!param](/NodalKernels/ClusterDynamicsNodalKernel/temperature)) is the temperature, $\sigma$ ([!param](/NodalKernels/ClusterDynamicsNodalKernel/sigma)) is the interfacial energy, $\Omega$ ([!param](/NodalKernels/ClusterDynamicsNodalKernel/Omega)) is the enthalpy term, and $\Delta S$ ([!param](/NodalKernels/ClusterDynamicsNodalKernel/DeltaS)) is the non-configurational entropy term.

### Intra-Variable Jacobian

The non-AD version (`ClusterDynamicsNodalKernel`) provides a full analytical Jacobian via intra-variable coupling using the `setJacobian` framework API.

+Monomer row ($i=0$):+

Diagonal entry (coupling of the monomer equation to $C_1$):

!equation
\frac{\partial F_0}{\partial C_1} = k_s + 4\beta_1 C_1 + \sum_{n \geq 2} \beta_n C_n

Off-diagonal entries (coupling of the monomer equation to $C_n$ for $n \geq 2$):

!equation
\frac{\partial F_0}{\partial C_n} = \beta_n C_1 - \mu_n \alpha_n, \quad n \geq 2, \quad \mu_n = \begin{cases} 2 & n = 2 \\ 1 & n > 2 \end{cases}

+Cluster rows (+$i \geq 1$, $n = i+1$+):+

Coupling to the monomer column:

!equation
\frac{\partial F_i}{\partial C_1} =
\begin{cases}
-2\beta_1 C_1 + \beta_2 C_2 & n = 2 \\
-\beta_{n-1} C_{n-1} + \beta_n C_n & n > 2
\end{cases}

This piecewise form is the derivative of the residual
$F_i = -(\dot{C}_n^{\text{in}} - \beta_n C_1 C_n + \alpha_{n+1} C_{n+1} - \alpha_n C_n)$
with respect to the monomer concentration $C_1$. The Jacobian expressions retain the same form for both rate models because $\beta_n$ and $\alpha_n$ depend on cluster size and user-specified parameters, but not on the solution variables themselves.

Lower-diagonal entry (coupling to cluster $n-1$, for $n > 2$):

!equation
\frac{\partial F_i}{\partial C_{n-1}} = -\beta_{n-1} C_1

Diagonal entry:

!equation
\frac{\partial F_i}{\partial C_n} = \beta_n C_1 + \alpha_n

Upper-diagonal entry (coupling to cluster $n+1$, for $n+1 \leq N$):

!equation
\frac{\partial F_i}{\partial C_{n+1}} = -\alpha_{n+1}

!alert note title=Problem Requirement
The intra-variable Jacobian capability requires `use_hash_table_matrix_assembly = true` see below input file [!ref](cd_inputfile)

## Example Input Syntax

!listing test/tests/cluster_dynamics/cluster_dynamics_50_combined.i block=Problem NodalKernels id=cd_inputfile

!syntax parameters /NodalKernels/ClusterDynamicsNodalKernel

!syntax inputs /NodalKernels/ClusterDynamicsNodalKernel

!syntax children /NodalKernels/ClusterDynamicsNodalKernel

## See Also

- [ClusterAverageRadius.md] - Average cluster radius postprocessor
- [ClusterTotalDensity.md] - Total cluster density postprocessor
- [ClusterSizeConcentration.md] - Single cluster size concentration postprocessor
