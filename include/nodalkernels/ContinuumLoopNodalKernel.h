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

#include "ArrayNodalKernel.h"
#include "ClusterContinuumFlux.h"
#include "ClusterContinuumLayout.h"

/**
 * Self-absorption and thermal emission of a family of two-dimensional loops whose size
 * distribution is discrete up to switch_size and a Fokker-Planck continuum above it (Pokor et al.
 * 2004, Section 4.1), the alternative to the L0/L1 grouping of
 * GenericInterstitialLoopGroupedNodalKernel. The rate coefficients are those of that kernel: loop
 * radius sqrt(n V_at / (pi b)), the size-dependent bias Z(n), and detailed-balance emission from a
 * capillary-law binding energy. The array variable has one component for the monomer, one for
 * each discrete size 2 .. switch_size, and one for each continuum cell (ClusterContinuumLayout).
 *
 * This version has no opposite-species terms, so it is the one-family model. Choose a largest size
 * N that leaves the distribution room to grow: the flux out of the last cell is zero, so a
 * population that reaches it piles up there.
 */
class ContinuumLoopNodalKernel : public ArrayNodalKernel
{
public:
  static InputParameters validParams();
  ContinuumLoopNodalKernel(const InputParameters & parameters);

protected:
  virtual void computeQpResidual(RealEigenVector & residual) override;
  virtual void computeQpJacobian() override;

  /// Loop absorption coefficient beta(n) for a monomer, 2 pi r(n) Z(n) D / V_at
  Real beta(unsigned int n) const;

  /// Emission coefficient alpha(n) = beta(n - 1) exp(-E_b(n) / kT), n >= 2
  Real alpha(unsigned int n) const;

  /// Fill the coefficients on first use
  void ensureCoefficients();

  const Real _generation;
  const Real _sink;
  const Real _temperature;
  const Real _atomic_volume;
  const Real _lattice_parameter;
  const Real _burgers_vector;
  const Real _diffusivity;
  const Real _z;
  const Real _z_loop;
  const Real _alpha_loop;
  const Real _formation_energy_eV;
  const Real _energy_exponent;
  const Real _energy_coefficient_eV;
  const Real _energy_constant_eV;

  const ClusterContinuumLayout _layout;
  const ClusterContinuumFlux _flux;

  ClusterContinuumCoefficients _coefficients;
  bool _have_coefficients = false;
  std::vector<Real> _rate;
  std::vector<ClusterContinuumFlux::Triplet> _triplets;
};
