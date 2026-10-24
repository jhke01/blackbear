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

#pragma once

#include "GeneralPostprocessor.h"

/**
 * Upper limit for the time step, for use as the timestep_limiting_postprocessor of a time stepper.
 * The limit keeps the relative change of a set of postprocessors (the resolvable cluster densities,
 * which rise by many orders of magnitude in a short interval of dose) below a target in one step:
 *
 *   limit = min(time_fraction * t + minimum_cap, dt * target / max(change, 1e-12)),
 *   change = max over the postprocessors of |v - v_old| / max(|v|, value_floor),
 *
 * where dt is the step just completed. At the initial time, and while nothing has changed, the limit
 * is the first term, which is also an upper bound on the step in proportion to the time.
 */
class RelativeChangeTimestepLimit : public GeneralPostprocessor
{
public:
  static InputParameters validParams();
  RelativeChangeTimestepLimit(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override {}
  virtual void finalize() override;
  virtual Real getValue() const override { return _value; }

protected:
  const Real _target;
  const Real _value_floor;
  const Real _time_fraction;
  const Real _minimum_cap;

  /// The monitored postprocessors at this step and at the previous one
  std::vector<const PostprocessorValue *> _current;
  std::vector<const PostprocessorValue *> _old;

  Real _value;
};
