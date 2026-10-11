# ADCascadeSpectrumNodalKernel / CascadeSpectrumNodalKernel

!syntax description /NodalKernels/CascadeSpectrumNodalKernel

## Description

See [theory/defect_dynamics.md] for a full derivation connecting this kernel's equation to the
rest of the defect-dynamics model.

`CascadeSpectrumNodalKernel` adds a constant generation source to one defect species' explicit
cluster sizes 2, 3, and 4 (array components 1, 2, 3 -- component 0 is the monomer), seeding
clusters directly from the cascade instead of letting the entire generation rate go to the
monomer row. It does not itself reduce the monomer's own generation term; the enclosing
[DefectDynamics](syntax/DefectDynamics/index.md) action does that on the core void/loop kernel's
own [!param](/NodalKernels/VacancyVoidGroupedNodalKernel/generation) parameter directly.

## Equation

Pokor et al. (2004) Eqs 1a-1e seed clusters of size 2, 3, and 4 directly from the cascade, with
each size's own fraction $f_n$ of the total generation rate $G$ (where $G$ is
[!param](/NodalKernels/CascadeSpectrumNodalKernel/generation),
$G = $ `dpa_rate * cascade_efficiency`). Applied literally, the paper's own
$G(n) = G f_n$ form counts $n$ point defects for every cluster seeded at size $n$, so the total
point-defect production rate becomes $1 + f_2 + 2f_3 + 3f_4$ times $G$ whenever any fraction is
nonzero -- more than the dpa rate actually produced. To keep total point-defect production
conserved at exactly $G$ regardless of how $f_2$, $f_3$ ([!param](/NodalKernels/CascadeSpectrumNodalKernel/f3)),
and $f_4$ ([!param](/NodalKernels/CascadeSpectrumNodalKernel/f4)) split it, each cluster size's
source is instead

!equation
\text{source}(n) = \frac{G f_n}{n}, \qquad n \in \{2, 3, 4\},

so that $n \times \text{source}(n) = G f_n$ contributes exactly its own fraction of the total
rate. This is a pure, state-independent, time-independent constant, so this kernel's Jacobian
contribution is exactly zero and it is AD-safe and provided in both AD and non-AD variants.

Cluster sizes 2, 3, and 4 are guaranteed to be explicit (ungrouped) array components whenever the
coupled array variable represents at least 4 cluster sizes, which this kernel's constructor
checks.

## Example Input Syntax

!listing test/tests/defect_dynamics/cascade_spectrum_small.i block=NodalKernels/vacancy_cascade

!syntax parameters /NodalKernels/CascadeSpectrumNodalKernel

!syntax inputs /NodalKernels/CascadeSpectrumNodalKernel

!syntax children /NodalKernels/CascadeSpectrumNodalKernel

## See Also

- [theory/defect_dynamics.md]
- [GenericVacancyVoidGroupedNodalKernel.md]
- [GenericInterstitialLoopGroupedNodalKernel.md]
- [GenericDislocationRecoveryNodalKernel.md]
