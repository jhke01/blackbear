/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*                       BlackBear                              */
/*                                                              */
/*           (c) 2017 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "RelativeChangeTimestepLimit.h"

registerMooseObject("BlackBearApp", RelativeChangeTimestepLimit);

InputParameters
RelativeChangeTimestepLimit::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addClassDescription(
      "Time-step limit that keeps the relative change of the given postprocessors below a target in "
      "one step, bounded above by time_fraction * t + minimum_cap. Use it as the "
      "timestep_limiting_postprocessor of a time stepper.");
  params.addRequiredParam<std::vector<PostprocessorName>>(
      "postprocessors", "Postprocessors whose relative change per step is limited.");
  params.addRangeCheckedParam<Real>(
      "relative_change_target",
      0.1,
      "relative_change_target > 0",
      "Largest relative change of a monitored postprocessor in one step. Smaller values follow a "
      "steep rise more closely and take more steps.");
  params.addRangeCheckedParam<Real>(
      "value_floor",
      1.0,
      "value_floor > 0",
      "A postprocessor below this value is measured against the floor, so that the noise of a "
      "population that has not yet formed does not limit the step.");
  params.addRangeCheckedParam<Real>("time_fraction",
                                    0.01,
                                    "time_fraction >= 0",
                                    "The step never exceeds time_fraction * t + minimum_cap.");
  params.addRangeCheckedParam<Real>(
      "minimum_cap", 0.5, "minimum_cap > 0", "See time_fraction. In the units of the time.");
  params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
  return params;
}

RelativeChangeTimestepLimit::RelativeChangeTimestepLimit(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _target(getParam<Real>("relative_change_target")),
    _value_floor(getParam<Real>("value_floor")),
    _time_fraction(getParam<Real>("time_fraction")),
    _minimum_cap(getParam<Real>("minimum_cap")),
    _value(0.0)
{
  const auto & names = getParam<std::vector<PostprocessorName>>("postprocessors");
  if (names.empty())
    paramError("postprocessors", "must name at least one postprocessor.");
  for (const auto i : index_range(names))
  {
    _current.push_back(&getPostprocessorValue("postprocessors", i));
    _old.push_back(&getPostprocessorValueOld("postprocessors", i));
  }
}

void
RelativeChangeTimestepLimit::finalize()
{
  const Real upper = _time_fraction * _t + _minimum_cap;

  // Old values of postprocessors are not stored at the initial step
  Real change = 0.0;
  if (_t_step > 0)
    for (const auto i : index_range(_current))
    {
      const Real value = *_current[i];
      change = std::max(change,
                        std::fabs(value - *_old[i]) / std::max(std::fabs(value), _value_floor));
    }

  _value = std::min(upper, _dt * _target / std::max(change, 1.0e-12));
}
