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

void
ClusterContinuumFlux::evaluate(const Real * state,
                               const ClusterContinuumCoefficients & coefficients,
                               std::vector<Flux> & flux) const
{
  const unsigned int sw = _layout.switchSize();
  const unsigned int nc = _layout.numCells();
  flux.assign(numFluxes(), Flux());

  const Real monomer = positive(state[0]);
  const Real monomer_derivative = positiveDerivative(state[0]);

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
  }

  // Slopes of the cell averages, limited so that no new extremum appears. The first and last
  // cell keep a zero slope.
  std::vector<Real> slope(nc, 0.0);
  if (_use_limiter)
    for (unsigned int k = 1; k + 1 < nc; ++k)
    {
      const Real previous = positive(state[_layout.cell(k - 1).component]);
      const Real current = positive(state[_layout.cell(k).component]);
      const Real next = positive(state[_layout.cell(k + 1).component]);
      const Real dx_left = _layout.cell(k).center - _layout.cell(k - 1).center;
      const Real dx_right = _layout.cell(k + 1).center - _layout.cell(k).center;
      const Real dx_center = _layout.cell(k + 1).center - _layout.cell(k - 1).center;
      slope[k] = minmod(limiter_theta * (current - previous) / dx_left,
                        (next - previous) / dx_center,
                        limiter_theta * (next - current) / dx_right);
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
  }
}

void
ClusterContinuumFlux::rate(const Real * state,
                           const ClusterContinuumCoefficients & coefficients,
                           Real generation,
                           Real sink,
                           std::vector<Real> & rate) const
{
  std::vector<Flux> flux;
  evaluate(state, coefficients, flux);

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
}

void
ClusterContinuumFlux::jacobian(const Real * state,
                               const ClusterContinuumCoefficients & coefficients,
                               Real sink,
                               std::vector<Triplet> & triplets) const
{
  std::vector<Flux> flux;
  evaluate(state, coefficients, flux);

  const unsigned int sw = _layout.switchSize();
  const unsigned int nc = _layout.numCells();
  triplets.clear();

  // Row r of -rate is minus the derivative of d C_r / dt; the monomer row is the same sum weighted
  // by the number of point defects each row holds, plus the sink.
  std::vector<Real> monomer_row(_layout.componentCount(), 0.0);
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
