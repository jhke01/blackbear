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

#include "ClusterContinuumLayout.h"

#include <algorithm>
#include <cmath>
#include <vector>

/**
 * Linear interpolation of a function of the integer size at a real-valued size, clamped to
 * 0 .. last, as in the SUNDIALS Pokor benchmark.
 */
template <typename F>
Real
interpolateClusterCoefficient(const F & f, Real size, int last)
{
  const int floor_size = static_cast<int>(std::floor(size));
  const Real frac = size - static_cast<Real>(floor_size);
  const int lo = std::min(std::max(floor_size, 0), last);
  const int hi = std::min(std::max(floor_size + 1, 0), last);
  return (1.0 - frac) * f(static_cast<unsigned int>(lo)) + frac * f(static_cast<unsigned int>(hi));
}

/**
 * Absorption and emission coefficients of one cluster family, evaluated where the continuum tail
 * needs them: at the discrete sizes, and at the centers and faces of the continuum cells. The
 * coefficients at a real-valued size are interpolated linearly between the two neighboring
 * integer sizes, as in the SUNDIALS Pokor benchmark.
 */
struct ClusterContinuumCoefficients
{
  /// beta(n), n = 1 .. switch_size (index 0 unused)
  std::vector<Real> beta_discrete;
  /// alpha(n), n = 2 .. switch_size + 1 (indices 0 and 1 unused)
  std::vector<Real> alpha_discrete;
  /// beta and alpha at the center of each cell
  std::vector<Real> beta_center, alpha_center;
  /// beta and alpha at the face between cell f - 1 and cell f, f = 1 .. numCells - 1 (index 0 unused)
  std::vector<Real> beta_face, alpha_face;

  /**
   * Fill the coefficients from functions of the integer size. beta(0) and alpha(0), alpha(1) must
   * return zero; beta(n) is the absorption rate coefficient of a cluster of n point defects for a
   * monomer, alpha(n) the rate at which a cluster of n loses one.
   */
  template <typename Beta, typename Alpha>
  void build(const ClusterContinuumLayout & layout, const Beta & beta, const Alpha & alpha)
  {
    const unsigned int sw = layout.switchSize();
    const unsigned int nc = layout.numCells();
    const int last = static_cast<int>(layout.numClusterSizes());

    auto interpolate = [last](const auto & f, Real size)
    { return interpolateClusterCoefficient(f, size, last); };

    beta_discrete.assign(sw + 1, 0.0);
    alpha_discrete.assign(sw + 2, 0.0);
    for (unsigned int n = 1; n <= sw; ++n)
      beta_discrete[n] = beta(n);
    for (unsigned int n = 2; n <= sw + 1; ++n)
      alpha_discrete[n] = alpha(n);

    beta_center.assign(nc, 0.0);
    alpha_center.assign(nc, 0.0);
    beta_face.assign(nc, 0.0);
    alpha_face.assign(nc, 0.0);
    for (unsigned int k = 0; k < nc; ++k)
    {
      beta_center[k] = interpolate(beta, layout.cell(k).center);
      alpha_center[k] = interpolate(alpha, layout.cell(k).center);
      if (k >= 1)
      {
        beta_face[k] = interpolate(beta, layout.cell(k - 1).right);
        alpha_face[k] = interpolate(alpha, layout.cell(k - 1).right);
      }
    }
  }
};

/**
 * Rate coefficients beta(n) at which clusters of one family absorb a monomer of the OTHER point
 * defect type, which shrinks them by one. Needed at the discrete sizes 1 .. switch_size + 1 and at
 * the centers and faces of the continuum cells.
 */
struct ClusterContinuumAbsorption
{
  /// beta(n), n = 1 .. switch_size + 1 (index 0 unused)
  std::vector<Real> discrete;
  /// beta at the center of each cell
  std::vector<Real> center;
  /// beta at the face between cell f - 1 and cell f, f = 1 .. numCells - 1 (index 0 unused)
  std::vector<Real> face;

  template <typename Beta>
  void build(const ClusterContinuumLayout & layout, const Beta & beta)
  {
    const unsigned int sw = layout.switchSize();
    const unsigned int nc = layout.numCells();
    const int last = static_cast<int>(layout.numClusterSizes());
    auto interpolate = [last](const auto & f, Real size)
    { return interpolateClusterCoefficient(f, size, last); };

    discrete.assign(sw + 2, 0.0);
    for (unsigned int n = 1; n <= sw + 1; ++n)
      discrete[n] = beta(n);
    center.assign(nc, 0.0);
    face.assign(nc, 0.0);
    for (unsigned int k = 0; k < nc; ++k)
    {
      center[k] = interpolate(beta, layout.cell(k).center);
      if (k >= 1)
        face[k] = interpolate(beta, layout.cell(k - 1).right);
    }
  }
};

/**
 * The other point-defect family that this one reacts with. Both families must use the same
 * ClusterContinuumLayout. Column numbers of the Jacobian at or above componentCount() belong to the
 * partner: column componentCount() + j is partner component j.
 */
struct ClusterContinuumOpposite
{
  /// This family absorbing the partner's monomer (shrinking by one)
  const ClusterContinuumAbsorption * own_absorption = nullptr;
  /// The partner's clusters absorbing this family's monomer
  const ClusterContinuumAbsorption * partner_absorption = nullptr;
  /// The partner's concentrations, componentCount() values
  const Real * partner_state = nullptr;
};

/**
 * Self-reaction fluxes and rates of one cluster family on a ClusterContinuumLayout, following
 * the discrete-plus-Fokker-Planck scheme of Pokor et al. (2004), Section 4.1, as implemented in
 * the SUNDIALS Pokor benchmark.
 *
 * Sizes 1 .. switch_size are individual unknowns. Above them the concentration per unit size is
 * averaged over finite-volume cells, and the flux between cells has an advective part (drift
 * beta*C1 - alpha, upwinded) and a diffusive part (half of beta*C1 + alpha). With use_limiter the
 * cell values at a face are reconstructed with a minmod-limited slope, as in SUNDIALS; without it
 * they are the cell averages, which is first-order upwind and has an exact Jacobian. The flux out
 * of the largest cell is zero. Concentrations enter through max(C, 0).
 *
 * Flux numbering: flux i < switch_size is the flux J_{i+1} from size i+1 to size i+2, so the last
 * discrete flux is the one into the first cell; flux switch_size - 1 + f is the flux across the face
 * between cells f - 1 and f. This class needs no MOOSE objects, so it can be compared directly
 * with the SUNDIALS right-hand side.
 */
class ClusterContinuumFlux
{
public:
  /// One entry of the Jacobian of the residual (minus the rate)
  struct Triplet
  {
    unsigned int row;
    unsigned int column;
    Real value;
  };

  ClusterContinuumFlux(const ClusterContinuumLayout & layout, bool use_limiter);

  unsigned int numFluxes() const { return _layout.switchSize() + _layout.numCells() - 1; }

  /**
   * d/dt of every component from the self reactions of this family only: the monomer gains
   * generation, loses sink*C1 and the point defects that the clusters take up.
   * @param state array of componentCount() concentrations
   */
  void rate(const Real * state,
            const ClusterContinuumCoefficients & coefficients,
            Real generation,
            Real sink,
            std::vector<Real> & rate,
            const ClusterContinuumOpposite * opposite = nullptr) const;

  /**
   * Jacobian of -rate with respect to the state. With use_limiter the derivative goes through the
   * minmod slopes: the slope of a cell depends on its two neighbors, and wherever minmod selects
   * one of its three arguments the derivative is that argument's, so the matrix is exact except
   * at the points where the selection changes. Entries of one row are consecutive and merged, and
   * row 0 (the monomer) is dense.
   */
  void jacobian(const Real * state,
                const ClusterContinuumCoefficients & coefficients,
                Real sink,
                std::vector<Triplet> & triplets,
                const ClusterContinuumOpposite * opposite = nullptr) const;

private:
  struct Flux
  {
    Real value = 0.0;
    unsigned int count = 0;

    /// The monomers (own and partner), the two cells next to the face and the cells their limited
    /// slopes depend on
    static constexpr unsigned int max_components = 8;
    unsigned int component[max_components] = {};
    Real partial[max_components] = {};

    void add(unsigned int c, Real d)
    {
      for (unsigned int i = 0; i < count; ++i)
        if (component[i] == c)
        {
          partial[i] += d;
          return;
        }
      component[count] = c;
      partial[count] = d;
      ++count;
    }
  };

  /// With opposite, the flux of the shrinking reactions is added to each flux
  void evaluate(const Real * state,
                const ClusterContinuumCoefficients & coefficients,
                const ClusterContinuumOpposite * opposite,
                std::vector<Flux> & flux) const;

  /// Sum over sizes 2 .. N of beta(n) C_n, the rate per unit partner monomer at which the clusters
  /// absorb it
  Real absorptionSum(const Real * state, const ClusterContinuumAbsorption & absorption) const;

  /// Index of the flux into cell k from the left (k = 0: the last discrete flux)
  unsigned int fluxIntoCell(unsigned int k) const { return _layout.switchSize() - 1 + k; }

  const ClusterContinuumLayout & _layout;
  const bool _use_limiter;
};
