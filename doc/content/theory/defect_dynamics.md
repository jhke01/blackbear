# Defect-Dynamics Theory: Vacancy Voids and Interstitial Loops

This page derives the governing equations behind the two-family defect-dynamics
model implemented by [GenericVacancyVoidGroupedNodalKernel.md],
[GenericInterstitialLoopGroupedNodalKernel.md], and
[PointDefectRecombinationNodalKernel.md]. It connects those three objects'
equations into one narrative so the physics is derived once rather than
repeated three times; the per-class pages remain the authoritative parameter
reference.

## Model Overview

Under irradiation, vacancies and self-interstitial atoms are produced in
pairs and either recombine, migrate to a fixed sink (e.g. dislocations), or
cluster into three-dimensional voids (vacancy clusters) or two-dimensional
loops (interstitial clusters) [!citep](pokor2004). This model tracks two
independent cluster-size distributions, $C_{v,n}$ (void concentration at
size $n$) and $C_{i,n}$ (loop concentration at size $n$), each using the
same explicit/grouped L0/L1 representation as the single-species
Cu-precipitation model; see [GenericClusterDynamicsGroupedNodalKernel.md]
for the grouping, moment evolution, and reconstruction-limiter derivation,
which both void and loop kernels reuse unchanged. That representation
keeps the first two moments ($L_0$, the bin-averaged concentration, and
$L_1$, its slope across the bin) of each grouped bin rather than one
representative concentration per bin, which is what lets the grouped
equations reproduce the ungrouped ones exactly in the limit of unit-width
bins and remain accurate for much wider bins than a single-moment grouping
would allow [!citep](golubov2001).

Each family has its own grouping (whether to group, how many explicit sizes, the bin spacing and
growth ratio), so a coarse grid for the voids and a fine grid for the loops is simply a different
value on each species; see [DefectDynamics/index.md].

Three physical processes couple or drive the two distributions:

- +Self-absorption and thermal emission+ -- a monomer of the same species is
  absorbed onto, or emitted from, an existing cluster. Implemented
  independently by [GenericVacancyVoidGroupedNodalKernel.md] (void family) and
  [GenericInterstitialLoopGroupedNodalKernel.md] (loop family).
- +Opposite-species absorption+ -- a monomer of the *opposite* species is
  absorbed onto an existing cluster, shrinking it by one unit. Implemented by
  the same two kernels, active only when each is given the other family's
  array variable.
- +Point-defect recombination+ -- a vacancy monomer and an interstitial
  monomer annihilate each other directly, before either one clusters.
  Implemented by [PointDefectRecombinationNodalKernel.md].

## Self-Absorption and Emission

For cluster size $n\geq 2$, the self-flux between sizes $n$ and $n+1$ follows
the same cluster-reaction form used throughout the cluster-dynamics kernels
in this code ([GenericClusterDynamicsNodalKernel.md]), and the same form
used for interstitial-loop and vacancy-cluster evolution in the cluster
dynamics literature [!citep](pokor2004):

!equation
J_n = \beta_n C_1 C_n - \alpha_{n+1} C_{n+1},

where $\beta_n$ is the absorption rate coefficient for a monomer joining a
size-$n$ cluster and $\alpha_n$ is the emission rate coefficient for a
size-$n$ cluster losing one unit. The two families differ in cluster
geometry, so $\beta_n$ takes a different functional form for each.

### Void (Vacancy) Absorption: Finite-Sink-Strength Correction

Voids are three-dimensional, so the void radius follows the same spherical
form as the Cu-precipitation model,

!equation
r_n = \left(\frac{3nV_{at}}{4\pi}\right)^{1/3},

with the vacancy-monomer diffusivity following an Arrhenius law,
$D_v = D_{0,v}\exp(-E_{m,v}/k_BT)$. Unlike the Cu kernel's uncorrected
$4\pi r_n D_v/V_{at}$ form, the void absorption coefficient includes a
finite-sink-strength correction,

!equation
\beta_n = \frac{4\pi r_n D_v}{\left(1+\dfrac{a_0}{r_n}\right)V_{at}},

where $a_0$ is the lattice parameter. This correction becomes significant
once $r_n$ approaches $a_0$ (i.e. for the smallest voids) and suppresses the
absorption rate relative to the uncorrected form.

### Loop (Interstitial) Absorption: Dislocation-Loop Bias

Loops are two-dimensional, so the loop radius follows

!equation
r_n = \sqrt{\frac{nV_{at}}{\pi b}},

where $b$ is the Burgers vector (defaulting to $a_0/\sqrt3$). In place of the
void's finite-sink-strength correction, loop absorption uses a size-dependent
dislocation-loop bias factor, reflecting a dislocation loop's enhanced
capture efficiency for interstitials relative to vacancies
[!citep](pokor2004):

!equation
\beta_n = 2\pi r_n Z_i(n) D_i / V_{at},
\qquad
Z_i(n) = Z_{i,0} + \left[\sqrt{\frac{b}{8\pi a_0}}Z_{loop,i} - Z_{i,0}\right]
n^{-\alpha_{loop,i}/2},

with the interstitial-monomer diffusivity $D_i = D_{0,i}\exp(-E_{m,i}/k_BT)$.
$Z_i(n)$ is $\sqrt{b/(8\pi a_0)}\,Z_{loop,i}$ at $n=1$ (about $0.152\,Z_{loop,i}$ for
$b=a_0/\sqrt3$) and approaches $Z_{i,0}$ as $n$ grows, so $Z_{i,0}$ is the
large-loop bias and $Z_{loop,i}$ sets the small-loop bias. The exponent
$\alpha_{loop,i}$ controls how fast the approach is. A constant bias $Z_i(n)=Z_{i,0}$
at every size therefore needs $Z_{loop,i}=Z_{i,0}/\sqrt{b/(8\pi a_0)}$.

### Capillary-Law Thermal Emission

Both families close the absorption/emission pair by detailed balance from a
capillary-law formation energy, parameterized identically in form for both
species ($x \in \{v,i\}$):

!equation
E_f(n) = A_x n^{p_x} + C_x,
\qquad
E_b(n) = E_{form,x} - \left[E_f(n)-E_f(n-1)\right],
\qquad
\alpha_n = \beta_{n-1}\exp\left(-\frac{E_b(n)}{k_BT}\right).

The void exponent is $p_v \approx 2/3$ (surface-area scaling for a 3D
cluster); the loop exponent $p_i$ is independently supplied since a 2D
cluster's edge energy need not scale the same way.

### Monomer and Grouped-Bin Equations

Written in terms of $J_n$, the cluster equations away from the monomer are
$dC_n/dt = J_{n-1}-J_n$, and the monomer equation is

!equation
\frac{dC_1}{dt} = G_1 - k_s C_1 - 2J_1 - \sum_{n=2}^{N-1} J_n,

matching the single-species form in [GenericClusterDynamicsNodalKernel.md]
exactly, with $G_1$ the monomer generation rate and $k_s$ the linear
dislocation-sink coefficient. The grouped-bin moment equations
($L_0$, $L_1$ per bin) are likewise unchanged from
[GenericClusterDynamicsGroupedNodalKernel.md]; both void and loop kernels
substitute their own $\beta_n$/$\alpha_n$ into that same closed-form
derivation.

## Opposite-Species Absorption

The processes above only exchange a cluster with monomers of its *own*
species. A second absorption channel lets a monomer of the *opposite*
species shrink an existing cluster by one unit (an interstitial absorbed by
a void shrinks the void; a vacancy absorbed by a loop shrinks the loop).
Without this term, nothing limits how large either population's clusters
can grow from the opposite species' absorption.

For a cluster of family $x$ at size $n\geq 2$, the opposite family $y$'s
monomer shrinks it to size $n-1$ at rate $\beta^{opp}_{x,n} C_{y,1} C_{x,n}$.
The two coefficients take the *same functional form* as that family's own
self-absorption $\beta_n$ above, but with the opposite species' diffusivity
(and, for loops, the opposite species' bias parameters) in place of the
self-absorption values:

!equation
\beta^{opp}_{v,n} = \frac{4\pi r_n D_i}{\left(1+\dfrac{a_0}{r_n}\right)V_{at}}
\qquad\text{(void absorbing an interstitial monomer)},

!equation
\beta^{opp}_{i,n} = 2\pi r_n Z_v(n) D_v / V_{at},
\qquad
Z_v(n) = Z_{v,0} + \left[\sqrt{\frac{b}{8\pi a_0}}Z_{loop,v} - Z_{v,0}\right]
n^{-\alpha_{loop,v}/2}
\qquad\text{(loop absorbing a vacancy monomer)}.

$Z_v(n)$ has the same functional form as $Z_i(n)$ above, but with independent
vacancy-side bias parameters ($Z_{v,0}$, $Z_{loop,v}$, $\alpha_{loop,v}$),
since there is no reason the loop's efficiency at capturing a vacancy should
match its efficiency at capturing an interstitial.

### Piece A: Own-Array Feedback

Define $K_x[n] = \beta^{opp}_{x,n}\,C_{y,1}\,C_{x,n}$, the opposite-shrink flux
out of size $n$ (zero for $n<2$: a monomer has no smaller size to shrink
into). Substituting the +combined flux+

!equation
\tilde J_x[n] = J_{x,n} - K_x[n+1]

for the plain self-flux $J_{x,n}$ everywhere it appears reproduces the
*exact same* closed-form residual used by the explicit-size and grouped-bin
equations above -- no other structural change is needed for those rows. Only
the monomer equation needs an explicit extra term, and that term is +not+
doubled the way self-dimerization is: a dimer-or-larger cluster shrinking to
a monomer adds back exactly one unit,

!equation
\frac{dC_{x,1}}{dt}\Big|^{\text{Piece A}} = +K_x[2].

### Piece B: Opposite-Monomer Drain

The flux absorbed in Piece A has to come from somewhere: the opposite
family's monomer pool. One absorption event removes one monomer from that pool and shrinks one
cluster of family $x$ by one unit, so the two must happen at the same rate. This loss is therefore
the sum, over family $x$'s +entire own array+, of the same $\beta^{opp}_{x,n}$ coefficients already
used in Piece A --

!equation
\left.\frac{dC_{y,1}}{dt}\right|^{\text{Piece B}} =
-C_{y,1}\sum_{n=2}^{N_x} \beta^{opp}_{x,n}\,C_{x,n}.

The kernel on variable $y$ does not own family $x$'s parameters, so it must be told family $x$'s
geometry. For a void family $x$ the coefficient is the void form above. For a loop family $x$ (for
example a vacancy modeled as a loop) it is the loop form $2\pi r_n Z(n) D_y/V_{at}$ with family $x$'s
own radius law and bias: [GenericInterstitialLoopGroupedNodalKernel.md] takes this from
`opposite_geometry = loop` together with the opposite family's Burgers vector and bias parameters.
Using the void form for a loop family drains the monomer pool at the wrong rate, and the total
number of defects is no longer conserved (a factor of about 2.8 to 7 in the rate at the sizes of
interest for Pokor's alloys).

Because a nodal kernel can only write to its own variable's residual row, the
kernel sitting on variable $y$ must recompute this sum itself, which is why
each kernel carries a second, fully independent `ClusterGroupingLayout`
describing the *coupled* family's grouping (`opposite_num_cluster_sizes`,
`opposite_explicit_max`, etc.) -- the two families' grids need not match, and
reconstructing $C_{x,n}$ from the coupled array requires its own layout.

### Why Non-AD Only

Piece B needs the coupled array's full contents, obtained through
`Coupleable::coupledArrayValue()`. That accessor returns a plain, non-AD
`ArrayVariableValue`, and `ADArrayNodalKernel::setJacobian()` unconditionally
errors -- so, exactly as for [PointDefectRecombinationNodalKernel.md], an AD
kernel could never actually assemble the resulting cross-variable Jacobian
block. Both [GenericVacancyVoidGroupedNodalKernel.md] and
[GenericInterstitialLoopGroupedNodalKernel.md] therefore hand-code the
opposite-absorption Jacobian (both the intra-variable part, folded into an
"effective $\alpha$" alongside the existing chain rule, and the genuinely
cross-variable part, a dense row against the coupled array) only in their
non-AD variant; their self-absorption terms remain AD-correct in both
variants.

When `opposite_coupled_var` is left unset, every term above evaluates to
exactly zero and both kernels reduce to the self-absorption-only model
described in the previous section.

## Point-Defect Recombination

Independently of clustering, a vacancy monomer and an interstitial monomer
can annihilate directly. [PointDefectRecombinationNodalKernel.md] adds this
loss to each family's monomer component at rate

!equation
R = k_{rec} C_{v,1} C_{i,1},
\qquad
k_{rec} = \frac{4\pi r_{rec}(D_i+D_v)}{V_{at}}.

Two symmetric instances are required: one subtracting $R$ from the vacancy
monomer equation (coupled to the interstitial variable), one subtracting $R$
from the interstitial monomer equation (coupled to the vacancy variable).
Like the opposite-absorption terms above, this is a genuinely cross-variable
term and is implemented non-AD only, for the same reason.

## Dislocation Recovery

Both families' monomer equations include a dislocation-sink loss term, $k_s C_1$, where $k_s =
D Z \rho$ (the core kernels' own `sink` parameter -- see
[GenericVacancyVoidGroupedNodalKernel.md]'s own documentation of its `sink` input). By default
$\rho$, the dislocation density, is a fixed constant. [GenericDislocationRecoveryNodalKernel.md]
replaces this with a time-evolving $\rho(t)$ for species where that module is active, reflecting
the experimental observation (Pokor et al. 2004) that the initial dislocation network of a
cold-worked alloy anneals out under irradiation via dislocation climb: a net flux of absorbed
point defects lets dislocations climb and annihilate each other, at a rate set by

!equation
\frac{d\rho}{dt} = -K b^2 \rho^{3/2},

whose closed-form solution,

!equation
\rho(t) = \left(\rho_0^{-1/2} + \frac{1}{2}K b^2 t\right)^{-2},

depends only on the current simulation time and three constants: the initial density $\rho_0$,
the recovery rate constant $K$, and the dislocation network's own Burgers vector $b$ -- a
distinct quantity from any loop species' own geometric Burgers vector. Neither $\rho$, $K$, nor
this $b$ carries a species subscript in Pokor's own equations, so one shared `[Recovery]`
sub-block (`DefectRecoveryAction`) supplies the same $\rho_0$/$K$/$b$ to a
[GenericDislocationRecoveryNodalKernel.md] instance on each species present, each using that
species' own $D$ and $Z$. Because $\rho(t)$ never depends on the solution itself, this term's
Jacobian is a plain diagonal entry, recomputed each call from the current simulation time rather
than read from a stored constant -- and because it needs no coupling to any other variable, the
kernel is AD-safe, unlike the opposite-absorption and recombination terms above.

!alert note
Pokor's own text quotes $K$ in units of m$^{-1}$ for the three reactor conditions it reports
(300, 970, 14400), but the governing equation above requires $K$ in m$^{-1}$s$^{-1}$ for
dimensional consistency. Treat $K$ as a fit parameter in self-consistent SI units, not a literal
copy of the paper's quoted numbers.

When no `[Recovery]` sub-block is given, a species' sink stays exactly the constant
`sink = D*Z*dislocation_density` term described above; the two are not combined -- Pokor's model
has exactly one sink term per species, so `[Recovery]`, when present, replaces that constant
contribution rather than adding to it.

## Cellular Network Dislocation Density id=cellular-network

With `model = cellular_network` in the `[Recovery]` sub-block ([DefectRecoveryAction.md]) the dislocation density is not a closed-form
function of time but follows the network model of [!cite](rest2006), which was derived for the cellular dislocation network that forms in
irradiated fuel. Interstitial loops climb to the network and add line length to it, and adjacent edge dislocations annihilate each other:

!equation
\frac{d\rho}{dt} = \pi v_l n_l - \frac{4 v_l}{d_g}\rho - v_l \frac{\sqrt{f(\nu)/\pi}}{C_A C_\rho}\rho^{3/2},
\qquad v_l = \frac{2 Z_{iv} D_i c_i}{b}, \quad f(\nu) = \frac{1 - \nu/2}{1 - \nu}.

Here $n_l$ is the interstitial loop density, $D_i c_i$ the interstitial diffusivity and monomer site fraction (the monomer of the loop species), $Z_{iv}$ the relative
interstitial/vacancy bias ([!param](/DefectDynamics/Recovery/DefectRecoveryAction/bias)), $b$ the Burgers vector of the network,
$d_g$ the grain size, $\nu$ the Poisson ratio, and $C_A C_\rho$ the geometric factor of the cells. The second term is the loss at grain boundaries (omitted when
[!param](/DefectDynamics/Recovery/DefectRecoveryAction/grain_size) is 0; the paper states it can be neglected for the conditions it considers), and the third term is the self-annihilation.
The dislocation sink of each species is $D Z \rho$ with this $\rho$, as in the previous section.

### Loop density at steady state

The loop density obeys $dn_l/dt = (\sqrt{2}/2) D_i c_i^2 / \Omega^{5/3} - v_l n_l / d_l$, where the first term is the nucleation of loops from di-interstitials and the second the loss of loops
to the network, with $\rho = \pi n_l d_l$ relating the loop diameter $d_l$ to the line length. Setting the left side to zero gives

!equation
n_l = \left[\frac{\sqrt{2}\, b\, c_i\, \rho}{4 \pi Z_{iv}\Omega^{5/3}}\right]^{1/2},

so $n_l \propto \sqrt{\rho}$. Substituting this into the equation for $\rho$ and writing $u = \sqrt{\rho}$, $n_l = k u$ with $k = [\sqrt{2} b c_i / (4 \pi Z_{iv}\Omega^{5/3})]^{1/2}$ gives

!equation
\frac{du}{dt} = \frac{\pi v_l k}{2} - \frac{2 v_l}{d_g} u - \frac{C v_l}{2} u^2, \qquad C = \frac{\sqrt{f(\nu)/\pi}}{C_A C_\rho},

a Riccati equation. For constant coefficients $du/dt = c - \beta u - a u^2 = -a(u - u_+)(u - u_-)$, and $w = (u - u_+)/(u - u_-)$ satisfies
$dw/dt = -\sqrt{\beta^2 + 4 a c}\, w$, so the solution over a step is exact. [CellularNetworkDislocationDensity.md] advances $\rho$ this way at the end of each step with $c_i$ at the mean of its values at the two ends
of the step. It cannot become negative. Because the density is read by the sink kernels as it stands at the end of the previous step, the coupling is lagged. For the cases tried,
changing [!param](/DefectDynamics/DefectDynamicsAction/timestep_relative_change) from 0.1 to 0.03 and 0.01 did not change the result by more than 0.01 percent in the density, so
the default is enough there; this was not tested where $c_i$ changes much faster.

### What has been checked, and what has not

All numbers are from one set of parameters (SA304 steel in the OSIRIS reactor, Pokor et al. Table 5 and 6), one machine, with $Z_{iv}$ and $C_A C_\rho$ assumed.

- The exact step agrees with an ODE integration at a relative tolerance of $10^{-12}$ to better than $3 \times 10^{-14}$.
- Against the same equations in the SUNDIALS Pokor benchmark with matched physics (no grain-boundary sink in the monomer balance, the 2/1 dimer bookkeeping), the density and the cluster
  observables agree within 0.2 percent over 0.01 to 10 dpa at 603 K. This checks the implementation of the equations, not their physics.
- At 603 K and 10 dpa, from an initial density of $10^{10}$ m$^{-2}$ with $Z_{iv} = 0.1$ and $C_A C_\rho = 2$, the density is $1.4 \times 10^{14}$ m$^{-2}$.
  Pokor's law gives $2 \times 10^{9}$ m$^{-2}$. The model changes the sign of the trend as well as the size.
- The values of $Z_{iv}$ (0.1 and 0.01 were used; [!cite](rest2006) gives $5 \times 10^{-4}$ for UO$_2$) and of $C_A C_\rho$ (2 and 3) are not known for steel. Results change by a factor of 10 or
  more between $Z_{iv} = 0.1$ and 0.01, and by a factor of 1.6 in the density at 10 dpa between $C_A C_\rho = 2$ and 3 (from an initial density of $10^{15}$ m$^{-2}$ at 603 K: $2.1 \times 10^{14}$ and $3.3 \times 10^{14}$ m$^{-2}$).
- At 673 K ($Z_{iv} = 0.1$) initial densities of $10^{10}$ and $10^{15}$ m$^{-2}$ both end at $3.6 \times 10^{14}$ and $3.7 \times 10^{14}$ m$^{-2}$ at 10 dpa.
  At 603 K the same two starts end at $1.4 \times 10^{14}$ and $2.1 \times 10^{14}$ m$^{-2}$ at 10 dpa, so they have not converged there. The model has not been run at 573 K or 623 K.
- No dislocation-density data of a steel were compared with the model.

### Line length and lattice strain from the loop population

[ClusterContinuumMoment.md] can report, for a loop family, the line length per volume $\sum_n 2\pi r_n C_n / V_{at}$ (the loop diameter in $\rho = \pi n_l d_l$) and
the relative change of the lattice parameter $\Delta a / a_0 = \pi b n_l d_l^2 / 12$ ([!cite](rest2006)). With $r_n^2 = n V_{at}/(\pi b)$ the second is $\tfrac{1}{3}\sum_n n C_n$, one third of the volume fraction of
the atoms in loops. It assumes that each loop atom adds its atomic volume to the crystal and ignores other sources of strain; the sign is the opposite for vacancy loops, for which it should not be used.

## Cascade Spectrum

Every monomer-only input to this model puts the entire generation rate $G$ (dpa_rate times
cascade_efficiency) into $G_1$, the monomer source term in the self-absorption section above.
Pokor et al. (2004) Eqs 1a-1e instead seed clusters of size 2, 3, and 4 directly from the cascade,
with each size's own fraction ($f_2$, $f_3$, $f_4$) of $G$ going to that size and the monomer
getting the remainder:

!equation
G_1 = G(1 - f_2 - f_3 - f_4), \qquad G(n) = G f_n, \quad n \in \{2, 3, 4\}.

Taken literally, $G(n) = G f_n$ seeds $n$ point defects at size $n$, so the total point-defect
production rate becomes $G_1 + 2 G(2) + 3 G(3) + 4 G(4) = G(1 + f_2 + 2f_3 + 3f_4)$ -- larger than
$G$ itself whenever any fraction is nonzero. [GenericCascadeSpectrumNodalKernel.md] instead
divides each size's source by its own size,

!equation
G(n) = \frac{G f_n}{n}, \qquad n \in \{2, 3, 4\},

so that $n \times G(n) = G f_n$ contributes exactly its own fraction of the total rate, keeping
total point-defect production conserved at exactly $G$ regardless of how the fractions split it.
Since $f_2$, $f_3$, and $f_4$ are genuinely per-species numbers in Pokor's own equations (distinct
values for vacancies and interstitials), they are set per species rather than shared; a species
activates cascade seeding simply by setting $f_2 + f_3 + f_4 > 0$ on its own sub-block, with no
separate switch. Like dislocation recovery above, this term depends only on constants, not on the
solution, so its Jacobian contribution is exactly zero and the kernel is AD-safe.

## Coupling to Solute Clustering

[SoluteClusterDynamics/index.md] describes a third population, solute clusters such as Cu precipitates,
handled by [GenericClusterDynamicsGroupedNodalKernel.md]. In a vacancy-mediated alloy the solute monomer
moves only by exchanging with a vacancy, so irradiation speeds its diffusion up in proportion to how
far the vacancy monomer concentration $C_{v,1}$ is above its thermal-equilibrium value
$C_v^{th} = 4\exp(-E_f/k_BT)$:

!equation
D_s = D_s^{0}\,f_{rad}\,\frac{C_{v,1}}{C_v^{th}},

where $D_s^{0}$ is the solute monomer diffusivity without irradiation and $f_{rad}$ the optional static
enhancement factor (leave it at 1 when the live ratio is used, or the effect is counted twice). The
ratio $C_{v,1}/C_v^{th}$ is the vacancy supersaturation. It is read from component 0 of the vacancy
species' array, so it does not depend on that species' grouping. The coupling is one-way: the vacancy
balance does not feel the solute. Because $D_s$ appears in every rate coefficient of the solute
equations, each solute equation depends on $C_{v,1}$; the kernel provides this off-diagonal Jacobian
block, which MOOSE uses only for variable pairs named in the preconditioner (see
[GenericClusterDynamicsNodalKernel.md]). The thermal concentration is evaluated at the solute block's
temperature, which must equal the temperature of the `[DefectDynamics]` block, and the coupling has an
effect only for `rate_model = interfacial_energy`, the model that uses the monomer diffusivity.

## Numerical Considerations id=numerical-considerations

Two numerical-accuracy pitfalls are easy to hit when setting up a new case
with this model.

### Cluster-Size Cap Must Leave Headroom

[!param](/NodalKernels/VacancyVoidGroupedNodalKernel/num_cluster_sizes)
truncates the tracked cluster-size space; both families' grouped bins must
fit within it. If a run's cluster population grows large enough that its
average size approaches this cap, probability mass piles up against the
truncation boundary and the solution becomes numerically unreliable. There
is no general formula for "large enough" -- it depends on dose, dose rate,
and temperature -- so a case should check its own final average cluster size
against its cap and increase the cap (and, if using geometric grouping, the
number of bins needed to reach it) whenever that average approaches a
sizable fraction of the cap.

### Ungrouped Species and the Cost of the Cross-Species Block

A species with `use_grouping = false` gives every cluster size its own unknown. The direct solve
and MOOSE's storage for the Jacobian grow much faster than linearly with the number of unknowns,
so the action refuses an ungrouped species with more than 3000 sizes. Declaring the coupling
between the two species' arrays in the preconditioner (needed for the exact cross-species
Jacobian of opposite-species absorption and of recombination) makes MOOSE allocate a block of
storage for every pair of components of the two arrays. In a measured case with the void and loop
arrays at 1000 sizes each (about 1070 components in all, with grouping), declaring the pair raised
the peak process memory from about 0.1 GB to about 2.0 GB, and an allocation profile put most of
the extra memory in MOOSE's `Assembly::init`, in vectors of dense matrices. Without
the cross-species block the Newton iteration uses an approximate Jacobian for the coupling
between the families. It still converges to the same residual, but takes more nonlinear iterations
in the cases tried. For production-size problems the choice is between the full Jacobian and the
memory.

### Timestep Resolution Through the Nucleation Transient

`IterationAdaptiveDT`'s default growth heuristic (grow `dt` whenever Newton
converges in fewer than `optimal_iterations`) has no awareness of how fast
the physical solution itself is changing. Because this model's capillary
nucleation rate turns on sharply while Newton continues to converge in only
a few iterations, a loosely tuned growth schedule can let `dt` grow much
faster than the nucleation transient can be resolved, producing too few
output points to capture it accurately even though every individual step
converged. A tighter growth schedule (smaller `optimal_iterations` and
`iteration_window`, a `growth_factor` closer to 1) is recommended whenever a
run must resolve a sharp nucleation transient, at the cost of more steps
(and correspondingly more wall time) to reach the same final dose.

## See Also

- [GenericVacancyVoidGroupedNodalKernel.md]
- [GenericInterstitialLoopGroupedNodalKernel.md]
- [PointDefectRecombinationNodalKernel.md]
- [GenericDislocationRecoveryNodalKernel.md]
- [GenericCascadeSpectrumNodalKernel.md]
- [GenericClusterDynamicsGroupedNodalKernel.md]
- [GenericClusterDynamicsNodalKernel.md]
- [SoluteClusterDynamics/index.md]
- [DefectDynamics/index.md]
- [GroupedArrayFunctionIndexIC.md]
