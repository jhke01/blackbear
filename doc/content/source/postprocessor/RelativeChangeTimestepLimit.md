# RelativeChangeTimestepLimit

!syntax description /Postprocessors/RelativeChangeTimestepLimit

## Description

`RelativeChangeTimestepLimit` computes an upper limit for the time step, to be used as the
`timestep_limiting_postprocessor` of a time stepper such as `IterationAdaptiveDT`. It keeps the relative change of the
postprocessors in [!param](/Postprocessors/RelativeChangeTimestepLimit/postprocessors) below
[!param](/Postprocessors/RelativeChangeTimestepLimit/relative_change_target) in one step:

\begin{equation}
\Delta t_{\max} = \min\left(f\, t + c, \; \Delta t \, \frac{\epsilon}{\max_i \left|v_i - v_i^{\mathrm{old}}\right| / \max(|v_i|, v_0)}\right),
\end{equation}

where $\Delta t$ is the step just completed, $\epsilon$ is the target, $f$ and $c$ are
[!param](/Postprocessors/RelativeChangeTimestepLimit/time_fraction) and
[!param](/Postprocessors/RelativeChangeTimestepLimit/minimum_cap), and $v_0$ is
[!param](/Postprocessors/RelativeChangeTimestepLimit/value_floor). The second term keeps a quantity that
rises by many orders of magnitude in a short interval from being followed with steps that are too large. The first term is a bound in proportion to the
time, so that the step still grows when nothing changes. At the initial time, where no old value exists, only the first term is used.

The cluster-dynamics actions add this object themselves, for the resolvable densities, when
`timestep_relative_change` is set (see [DefectDynamics](syntax/DefectDynamics/index.md) and
[SoluteClusterDynamics](syntax/SoluteClusterDynamics/index.md)).

## Example Input Syntax

In this example the time step of a cluster-dynamics problem is limited by the relative change of the resolvable densities.

!listing test/tests/defect_dynamics/defect_dynamics_timestep_limit_action.i block=DefectDynamics

!syntax parameters /Postprocessors/RelativeChangeTimestepLimit

!syntax inputs /Postprocessors/RelativeChangeTimestepLimit

!syntax children /Postprocessors/RelativeChangeTimestepLimit
