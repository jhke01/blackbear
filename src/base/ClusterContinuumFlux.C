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

#include "ClusterContinuumFlux.h"

#include <array>

namespace
{
/// Concentrations enter the fluxes as max(C, 0), as in the SUNDIALS Pokor benchmark
inline Real
positive(Real value)
{
  return value > 0.0 ? value : 0.0;
}

/// Derivative of positive()
inline Real
positiveDerivative(Real value)
{
  return value > 0.0 ? 1.0 : 0.0;
}

/// Parameter of the generalized minmod limiter: 1 is minmod, 2 is the monotonized central limiter.
/// 1.5 is the value used by the SUNDIALS Pokor benchmark.
constexpr Real limiter_theta = 1.5;

/// Generalized minmod of the three one-sided slope estimates
inline Real
minmod(Real a, Real b, Real c)
{
  if (a > 0.0 && b > 0.0 && c > 0.0)
    return std::min(a, std::min(b, c));
  if (a < 0.0 && b < 0.0 && c < 0.0)
    return std::max(a, std::max(b, c));
  return 0.0;
}
} // namespace

ClusterContinuumFlux::ClusterContinuumFlux(const ClusterContinuumLayout & layout, bool use_limiter)
  : _layout(layout), _use_limiter(use_limiter)
{
}

Real
ClusterContinuumFlux::absorptionSum(const Real * state,
                                    const ClusterContinuumAbsorption & absorption) const
{
  Real sum = 0.0;
  for (unsigned int n = 2; n <= _layout.switchSize(); ++n)
    sum += absorption.discrete[n] * positive(state[n - 1]);
  for (unsigned int k = 0; k < _layout.numCells(); ++k)
    sum += absorption.center[k] * positive(state[_layout.cell(k).component]) * _layout.cell(k).width;
  return sum;
}

void
ClusterContinuumFlux::evaluate(const Real * state,
                               const ClusterContinuumCoefficients & coefficients,
                               const ClusterContinuumOpposite * opposite,
                               std::vector<Flux> & flux) const
{
  const unsigned int sw = _layout.switchSize();
  const unsigned int nc = _layout.numCells();
  flux.assign(numFluxes(), Flux());

  const Real monomer = positive(state[0]);
  const Real monomer_derivative = positiveDerivative(state[0]);

  // The partner's monomer and the column of its Jacobian entries
  const unsigned int partner_column = _layout.componentCount();
  const Real partner_monomer = opposite ? positive(opposite->partner_state[0]) : 0.0;
  const Real partner_monomer_derivative =
      opposite ? positiveDerivative(opposite->partner_state[0]) : 0.0;

  // J_n = beta(n) C_1 C_n - alpha(n + 1) C_{n + 1} between the discrete sizes n and n + 1. The
  // size n + 1 sits in component n, which is the first cell for n = switch_size.
  for (unsigned int n = 1; n <= sw; ++n)
  {
    Flux & f = flux[n - 1];
    const Real beta = coefficients.beta_discrete[n];
    const Real alpha = coefficients.alpha_discrete[n + 1];
    const Real c_n = positive(state[n - 1]);
    const Real c_next = positive(state[n]);
    f.value = beta * monomer * c_n - alpha * c_next;
    f.add(0, beta * c_n * monomer_derivative);
    f.add(n - 1, beta * monomer * positiveDerivative(state[n - 1]));
    f.add(n, -alpha * positiveDerivative(state[n]));

    // A cluster of size n + 1 that absorbs a partner monomer becomes size n: flux
    // -beta_opposite(n + 1) C_partner,1 C_{n + 1} from n + 1 down to n
    if (opposite)
    {
      const Real beta_opposite = opposite->own_absorption->discrete[n + 1];
      f.value -= beta_opposite * partner_monomer * c_next;
      f.add(n, -beta_opposite * partner_monomer * positiveDerivative(state[n]));
      f.add(partner_column, -beta_opposite * c_next * partner_monomer_derivative);
    }
  }

  // Slopes of the cell averages, limited so that no new extremum appears, and their derivatives
  // with respect to the cell averages of cells k - 1, k and k + 1. The first and last cell keep a
  // zero slope.
  std::vector<Real> slope(nc, 0.0);
  std::vector<std::array<Real, 3>> slope_partial(nc, std::array<Real, 3>{0.0, 0.0, 0.0});
  if (_use_limiter)
    for (unsigned int k = 1; k + 1 < nc; ++k)
    {
      const Real previous = positive(state[_layout.cell(k - 1).component]);
      const Real current = positive(state[_layout.cell(k).component]);
      const Real next = positive(state[_layout.cell(k + 1).component]);
      const Real dx_left = _layout.cell(k).center - _layout.cell(k - 1).center;
      const Real dx_right = _layout.cell(k + 1).center - _layout.cell(k).center;
      const Real dx_center = _layout.cell(k + 1).center - _layout.cell(k - 1).center;

      // The three one-sided estimates and their derivatives with respect to previous, current, next
      const std::array<Real, 3> estimate = {limiter_theta * (current - previous) / dx_left,
                                            (next - previous) / dx_center,
                                            limiter_theta * (next - current) / dx_right};
      const std::array<std::array<Real, 3>, 3> estimate_partial = {
          std::array<Real, 3>{-limiter_theta / dx_left, limiter_theta / dx_left, 0.0},
          std::array<Real, 3>{-1.0 / dx_center, 0.0, 1.0 / dx_center},
          std::array<Real, 3>{0.0, -limiter_theta / dx_right, limiter_theta / dx_right}};

      slope[k] = minmod(estimate[0], estimate[1], estimate[2]);
      if (slope[k] != 0.0)
      {
        // minmod returns the estimate of smallest magnitude when all three agree in sign
        unsigned int selected = 0;
        for (unsigned int i = 1; i < 3; ++i)
          if (std::fabs(estimate[i]) < std::fabs(estimate[selected]))
            selected = i;
        const Real chain[3] = {positiveDerivative(state[_layout.cell(k - 1).component]),
                               positiveDerivative(state[_layout.cell(k).component]),
                               positiveDerivative(state[_layout.cell(k + 1).component])};
        for (unsigned int i = 0; i < 3; ++i)
          slope_partial[k][i] = estimate_partial[selected][i] * chain[i];
      }
    }

  // Flux across the face between cell f - 1 and cell f: the drift beta C_1 - alpha carries the
  // upwind cell value at the face, and the diffusion (beta C_1 + alpha) / 2 acts on the cell
  // averages. The slopes are held fixed in the derivatives.
  for (unsigned int face = 1; face < nc; ++face)
  {
    const ClusterContinuumCell & left = _layout.cell(face - 1);
    const ClusterContinuumCell & right = _layout.cell(face);
    const Real c_left = state[left.component];
    const Real c_right = state[right.component];
    const Real x_face = left.right;

    const Real left_inner = positive(c_left) + slope[face - 1] * (x_face - left.center);
    const Real right_inner = positive(c_right) + slope[face] * (x_face - right.center);
    const Real left_value = positive(left_inner);
    const Real right_value = positive(right_inner);
    const Real left_value_derivative = positiveDerivative(left_inner) * positiveDerivative(c_left);
    const Real right_value_derivative =
        positiveDerivative(right_inner) * positiveDerivative(c_right);

    const Real drift = coefficients.beta_face[face] * monomer - coefficients.alpha_face[face];
    const bool from_left = drift >= 0.0;
    const Real advective = from_left ? drift * left_value : drift * right_value;

    const Real spacing = right.center - left.center;
    const Real twice_diffusion_left =
        coefficients.beta_center[face - 1] * monomer + coefficients.alpha_center[face - 1];
    const Real twice_diffusion_right =
        coefficients.beta_center[face] * monomer + coefficients.alpha_center[face];
    const Real diffusive =
        -0.5 * (twice_diffusion_right * positive(c_right) - twice_diffusion_left * positive(c_left)) /
        spacing;

    Flux & f = flux[fluxIntoCell(face)];
    f.value = advective + diffusive;

    // The upwind face value depends on the cell it comes from and on the slope of that cell, which
    // depends on that cell and its two neighbors
    {
      const unsigned int upwind = from_left ? face - 1 : face;
      const Real inner_derivative =
          from_left ? positiveDerivative(left_inner) : positiveDerivative(right_inner);
      const Real offset = x_face - (from_left ? left.center : right.center);
      for (unsigned int i = 0; i < 3; ++i)
      {
        const int neighbor = static_cast<int>(upwind) - 1 + static_cast<int>(i);
        if (neighbor >= 0 && neighbor < static_cast<int>(nc) && slope_partial[upwind][i] != 0.0)
          f.add(_layout.cell(neighbor).component,
                drift * inner_derivative * offset * slope_partial[upwind][i]);
      }
    }

    f.add(0,
          monomer_derivative *
              (coefficients.beta_face[face] * (from_left ? left_value : right_value) -
               0.5 *
                   (coefficients.beta_center[face] * positive(c_right) -
                    coefficients.beta_center[face - 1] * positive(c_left)) /
                   spacing));
    f.add(left.component,
          (from_left ? drift * left_value_derivative : 0.0) +
              0.5 * twice_diffusion_left * positiveDerivative(c_left) / spacing);
    f.add(right.component,
          (from_left ? 0.0 : drift * right_value_derivative) -
              0.5 * twice_diffusion_right * positiveDerivative(c_right) / spacing);

    // Shrinking by the partner's monomers: drift -beta_opposite C_partner,1 (always toward smaller
    // sizes, so the face value comes from the right cell) and diffusion beta_opposite C_partner,1 / 2.
    // It has its own upwinding and shares the slopes.
    if (opposite)
    {
      const ClusterContinuumAbsorption & absorption = *opposite->own_absorption;
      const Real drift_opposite = -absorption.face[face] * partner_monomer;
      const Real twice_diffusion_left_opposite = absorption.center[face - 1] * partner_monomer;
      const Real twice_diffusion_right_opposite = absorption.center[face] * partner_monomer;

      f.value += drift_opposite * right_value -
                 0.5 *
                     (twice_diffusion_right_opposite * positive(c_right) -
                      twice_diffusion_left_opposite * positive(c_left)) /
                     spacing;

      const Real offset = x_face - right.center;
      for (unsigned int i = 0; i < 3; ++i)
      {
        const int neighbor = static_cast<int>(face) - 1 + static_cast<int>(i);
        if (neighbor >= 0 && neighbor < static_cast<int>(nc) && slope_partial[face][i] != 0.0)
          f.add(_layout.cell(neighbor).component,
                drift_opposite * positiveDerivative(right_inner) * offset * slope_partial[face][i]);
      }
      f.add(left.component,
            0.5 * twice_diffusion_left_opposite * positiveDerivative(c_left) / spacing);
      f.add(right.component,
            drift_opposite * right_value_derivative -
                0.5 * twice_diffusion_right_opposite * positiveDerivative(c_right) / spacing);
      f.add(partner_column,
            partner_monomer_derivative *
                (-absorption.face[face] * right_value -
                 0.5 *
                     (absorption.center[face] * positive(c_right) -
                      absorption.center[face - 1] * positive(c_left)) /
                     spacing));
    }
  }
}

void
ClusterContinuumFlux::rate(const Real * state,
                           const ClusterContinuumCoefficients & coefficients,
                           Real generation,
                           Real sink,
                           std::vector<Real> & rate,
                           const ClusterContinuumOpposite * opposite) const
{
  std::vector<Flux> flux;
  evaluate(state, coefficients, opposite, flux);

  const unsigned int sw = _layout.switchSize();
  const unsigned int nc = _layout.numCells();
  rate.assign(_layout.componentCount(), 0.0);

  // The point defects that the clusters take up are the monomer's loss: sum of n dC_n/dt over the
  // discrete sizes from 2 up, and of the first moment of the cell over the cells. A monomer pair
  // joining into a dimer therefore costs two monomers, as in the 2/1 bookkeeping of the grouped
  // kernels.
  Real defect_rate = 0.0;
  for (unsigned int n = 2; n <= sw; ++n)
  {
    const Real r = flux[n - 2].value - flux[n - 1].value;
    rate[n - 1] = r;
    defect_rate += static_cast<Real>(n) * r;
  }
  for (unsigned int k = 0; k < nc; ++k)
  {
    const Real into = flux[fluxIntoCell(k)].value;
    const Real out = k + 1 < nc ? flux[fluxIntoCell(k + 1)].value : 0.0;
    const Real r = (into - out) / _layout.cell(k).width;
    rate[_layout.cell(k).component] = r;
    defect_rate += _layout.cell(k).first_moment * r;
  }
  rate[0] = generation - sink * positive(state[0]) - defect_rate;

  if (opposite)
  {
    // The clusters that shrink by absorbing a partner monomer annihilate it, which costs the
    // partner's monomer a point defect each, and the partner's clusters absorbing this family's
    // monomer take that monomer
    rate[0] -= positive(opposite->partner_state[0]) * absorptionSum(state, *opposite->own_absorption) +
               positive(state[0]) * absorptionSum(opposite->partner_state, *opposite->partner_absorption);
  }
}

void
ClusterContinuumFlux::jacobian(const Real * state,
                               const ClusterContinuumCoefficients & coefficients,
                               Real sink,
                               std::vector<Triplet> & triplets,
                               const ClusterContinuumOpposite * opposite) const
{
  std::vector<Flux> flux;
  evaluate(state, coefficients, opposite, flux);

  const unsigned int sw = _layout.switchSize();
  const unsigned int nc = _layout.numCells();
  triplets.clear();

  // Row r of -rate is minus the derivative of d C_r / dt; the monomer row is the same sum weighted
  // by the number of point defects each row holds, plus the sink.
  const unsigned int partner_column = _layout.componentCount();
  std::vector<Real> monomer_row(opposite ? 2 * partner_column : partner_column, 0.0);
  auto addRow = [&](unsigned int row, Real defects, Real scale, const Flux & into, const Flux * out)
  {
    for (unsigned int i = 0; i < into.count; ++i)
    {
      triplets.push_back({row, into.component[i], -scale * into.partial[i]});
      monomer_row[into.component[i]] += defects * scale * into.partial[i];
    }
    if (out)
      for (unsigned int i = 0; i < out->count; ++i)
      {
        triplets.push_back({row, out->component[i], scale * out->partial[i]});
        monomer_row[out->component[i]] -= defects * scale * out->partial[i];
      }
  };

  for (unsigned int n = 2; n <= sw; ++n)
    addRow(n - 1, static_cast<Real>(n), 1.0, flux[n - 2], &flux[n - 1]);
  for (unsigned int k = 0; k < nc; ++k)
    addRow(_layout.cell(k).component,
           _layout.cell(k).first_moment,
           1.0 / _layout.cell(k).width,
           flux[fluxIntoCell(k)],
           k + 1 < nc ? &flux[fluxIntoCell(k + 1)] : nullptr);

  monomer_row[0] += sink * positiveDerivative(state[0]);

  if (opposite)
  {
    // Derivative of the monomer's loss to shrinking and to the partner's absorption (see rate())
    const Real own_monomer = positive(state[0]);
    const Real partner_monomer = positive(opposite->partner_state[0]);
    const ClusterContinuumAbsorption & own = *opposite->own_absorption;
    const ClusterContinuumAbsorption & partner = *opposite->partner_absorption;
    const Real * partner_state = opposite->partner_state;

    monomer_row[0] += positiveDerivative(state[0]) * absorptionSum(partner_state, partner);
    monomer_row[partner_column] +=
        positiveDerivative(partner_state[0]) * absorptionSum(state, own);
    for (unsigned int n = 2; n <= sw; ++n)
    {
      monomer_row[n - 1] += partner_monomer * own.discrete[n] * positiveDerivative(state[n - 1]);
      monomer_row[partner_column + n - 1] +=
          own_monomer * partner.discrete[n] * positiveDerivative(partner_state[n - 1]);
    }
    for (unsigned int k = 0; k < nc; ++k)
    {
      const unsigned int c = _layout.cell(k).component;
      const Real width = _layout.cell(k).width;
      monomer_row[c] += partner_monomer * own.center[k] * width * positiveDerivative(state[c]);
      monomer_row[partner_column + c] +=
          own_monomer * partner.center[k] * width * positiveDerivative(partner_state[c]);
    }
  }
  for (unsigned int c = 0; c < monomer_row.size(); ++c)
    if (monomer_row[c] != 0.0 || c == 0)
      triplets.push_back({0, c, monomer_row[c]});

  // Order by row and column and merge entries that fall on the same position
  std::sort(triplets.begin(),
            triplets.end(),
            [](const Triplet & a, const Triplet & b)
            { return a.row != b.row ? a.row < b.row : a.column < b.column; });
  std::vector<Triplet> merged;
  for (const auto & t : triplets)
    if (!merged.empty() && merged.back().row == t.row && merged.back().column == t.column)
      merged.back().value += t.value;
    else
      merged.push_back(t);
  triplets.swap(merged);
}
