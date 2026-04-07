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

#include "GenericClusterDynamicsNodalKernel.h"
#include "FEProblemBase.h"

#include "libmesh/libmesh_common.h"

namespace
{
constexpr Real kB = 1.380649e-23;
}

registerMooseObject("BlackBearApp", ClusterDynamicsNodalKernel);
registerMooseObject("BlackBearApp", ADClusterDynamicsNodalKernel);

template <bool is_ad>
InputParameters
GenericClusterDynamicsNodalKernelTempl<is_ad>::validParams()
{
  InputParameters params = GenericArrayNodalKernel<is_ad>::validParams();
  MooseEnum rate_model("simple interfacial_energy", "simple");
  MooseEnum diffusivity_model("constant arrhenius", "constant");
  params.addClassDescription(
      "Cluster dynamics nodal kernel for all cluster sizes (1 through N) in a single "
      "array variable. Array index i corresponds to cluster size n = i+1.");
  params.addRequiredParam<Real>("generation", "Monomer generation rate G_1");
  params.addRequiredParam<Real>("sink", "Linear sink coefficient k_s for monomer loss");
  params.addParam<MooseEnum>(
      "rate_model",
      rate_model,
      "Rate-coefficient model. 'simple' uses beta_n = beta0*n^(1/3) and "
      "alpha_n = alpha0*n^(1/3); 'interfacial_energy' derives coefficients from "
      "cluster geometry, interfacial energy, and detailed balance.");
  params.addParam<MooseEnum>(
      "diffusivity_model",
      diffusivity_model,
      "Diffusivity model used only when rate_model = interfacial_energy. "
      "'constant' uses monomer_diffusivity directly, while 'arrhenius' computes "
      "the diffusivity from D0, Q, and temperature.");
  params.addParam<Real>(
      "beta0", 0.0, "Base absorption coefficient for rate_model = simple.");
  params.addParam<Real>(
      "alpha0", 0.0, "Base emission coefficient for rate_model = simple.");
  params.addParam<Real>(
      "temperature", 0.0, "Temperature T [K] for the interfacial-energy rate model.");
  params.addParam<Real>(
      "monomer_diffusivity",
      0.0,
      "Monomer diffusion coefficient [m^2/s] for diffusivity_model = constant.");
  params.addParam<Real>(
      "D0", 0.0, "Diffusion prefactor D0 [m^2/s] for diffusivity_model = arrhenius.");
  params.addParam<Real>(
      "Q", 0.0, "Activation energy Q [J] for diffusivity_model = arrhenius.");
  params.addParam<Real>(
      "sigma", 0.0, "Interfacial energy sigma [J/m^2] for the interfacial-energy model.");
  params.addParam<Real>(
      "atomic_volume", 0.0, "Atomic volume V_at [m^3] for the interfacial-energy model.");
  params.addParam<Real>(
      "Omega", 0.0, "Enthalpy term Omega [J] for the interfacial-energy model.");
  params.addParam<Real>(
      "DeltaS",
      0.0,
      "Non-configurational entropy term DeltaS [J/K] for the interfacial-energy model.");
  return params;
}

template <bool is_ad>
GenericClusterDynamicsNodalKernelTempl<is_ad>::GenericClusterDynamicsNodalKernelTempl(
    const InputParameters & parameters)
  : GenericArrayNodalKernel<is_ad>(parameters),
    _generation(this->template getParam<Real>("generation")),
    _sink(this->template getParam<Real>("sink")),
    _rate_model(this->template getParam<MooseEnum>("rate_model")
                    .template getEnum<RateModel>()),
    _beta0(this->template getParam<Real>("beta0")),
    _alpha0(this->template getParam<Real>("alpha0")),
    _temperature(this->template getParam<Real>("temperature")),
    _monomer_diffusivity(this->template getParam<Real>("monomer_diffusivity")),
    _diffusivity_model(this->template getParam<MooseEnum>("diffusivity_model")
                           .template getEnum<DiffusivityModel>()),
    _D0(this->template getParam<Real>("D0")),
    _Q(this->template getParam<Real>("Q")),
    _sigma(this->template getParam<Real>("sigma")),
    _atomic_volume(this->template getParam<Real>("atomic_volume")),
    _Omega(this->template getParam<Real>("Omega")),
    _DeltaS(this->template getParam<Real>("DeltaS"))
{
  if (!(this->_fe_problem.useHashTableMatrixAssembly()))
    mooseError("ClusterDynamicsNodalKernel requires Problem/use_hash_table_matrix_assembly = true");

  if (_rate_model == RateModel::SIMPLE)
  {
    if (_beta0 <= 0.0)
      mooseError("ClusterDynamicsNodalKernel with rate_model = simple requires beta0 > 0.");
    if (_alpha0 <= 0.0)
      mooseError("ClusterDynamicsNodalKernel with rate_model = simple requires alpha0 > 0.");
  }
  else
  {
    if (_temperature <= 0.0)
      mooseError(
          "ClusterDynamicsNodalKernel with rate_model = interfacial_energy requires temperature > 0.");
    if (_diffusivity_model == DiffusivityModel::CONSTANT)
    {
      if (_monomer_diffusivity <= 0.0)
        mooseError("ClusterDynamicsNodalKernel with rate_model = interfacial_energy and "
                   "diffusivity_model = constant requires monomer_diffusivity > 0.");
    }
    else
    {
      if (_D0 <= 0.0)
        mooseError("ClusterDynamicsNodalKernel with rate_model = interfacial_energy and "
                   "diffusivity_model = arrhenius requires D0 > 0.");
      if (_Q <= 0.0)
        mooseError("ClusterDynamicsNodalKernel with rate_model = interfacial_energy and "
                   "diffusivity_model = arrhenius requires Q > 0.");
    }
    if (_sigma <= 0.0)
      mooseError(
          "ClusterDynamicsNodalKernel with rate_model = interfacial_energy requires sigma > 0.");
    if (_atomic_volume <= 0.0)
      mooseError(
          "ClusterDynamicsNodalKernel with rate_model = interfacial_energy requires atomic_volume > 0.");
  }
}

template <bool is_ad>
Real
GenericClusterDynamicsNodalKernelTempl<is_ad>::beta(const unsigned int n) const
{
  if (_rate_model == RateModel::SIMPLE)
    return _beta0 * std::cbrt(static_cast<Real>(n));

  return 4.0 * libMesh::pi * (radius(1) + radius(n)) * monomerDiffusivity() / atomicVolume();
}

template <bool is_ad>
Real
GenericClusterDynamicsNodalKernelTempl<is_ad>::alpha(const unsigned int n) const
{
  if (_rate_model == RateModel::SIMPLE)
    return _alpha0 * std::cbrt(static_cast<Real>(n));

  if (n <= 1)
    mooseError("alpha(n) is only defined for cluster sizes n >= 2 in interfacial_energy mode.");

  return beta(n - 1) * std::exp(-bindingEnergy(n) / (kB * _temperature));
}

template <bool is_ad>
Real
GenericClusterDynamicsNodalKernelTempl<is_ad>::atomicVolume() const
{
  return _atomic_volume;
}

template <bool is_ad>
Real
GenericClusterDynamicsNodalKernelTempl<is_ad>::monomerDiffusivity() const
{
  if (_diffusivity_model == DiffusivityModel::CONSTANT)
    return _monomer_diffusivity;

  return _D0 * std::exp(-_Q / (kB * _temperature));
}

template <bool is_ad>
Real
GenericClusterDynamicsNodalKernelTempl<is_ad>::radius(const unsigned int n) const
{
  return std::cbrt(3.0 * static_cast<Real>(n) * atomicVolume() / (4.0 * libMesh::pi));
}

template <bool is_ad>
Real
GenericClusterDynamicsNodalKernelTempl<is_ad>::bindingEnergy(const unsigned int n) const
{
  if (n <= 1)
    mooseError(
        "bindingEnergy(n) is only defined for cluster sizes n >= 2 in interfacial_energy mode.");

  const Real surface_prefactor =
      std::cbrt(36.0 * libMesh::pi) * std::pow(atomicVolume(), 2.0 / 3.0) * _sigma;
  return _Omega - _temperature * _DeltaS -
         surface_prefactor *
             (std::pow(static_cast<Real>(n), 2.0 / 3.0) -
              std::pow(static_cast<Real>(n - 1), 2.0 / 3.0));
}

template <bool is_ad>
void
GenericClusterDynamicsNodalKernelTempl<is_ad>::computeQpResidual(
    GenericRealEigenVector<is_ad> & residual)
{
  const auto n_comp = _u[_qp].size();
  residual.resize(n_comp);

  // c(i) = concentration of cluster of size n = i+1
  const auto c1 = _u[_qp](0); // monomer concentration

  // Component 0: monomer rate equation
  // dC_1/dt = G_1 - k_s*C_1
  //         - 2*beta(1)*C_1^2                        [monomer+monomer -> dimer]
  //         - sum_{i=1}^{N-1} beta(i+1)*C_1*C_{i+1} [monomer absorbed by larger cluster]
  //         + 2*alpha(2)*C_2                          [dimer dissociation]
  //         + sum_{i=2}^{N-1} alpha(i+1)*C_{i+1}     [larger cluster emission]
  {
    auto absorption = 2.0 * beta(1) * c1 * c1;
    GenericReal<is_ad> emission = 0.0;
    for (auto j = 1; j < n_comp; ++j)
    {
      const unsigned int nj = j + 1; // cluster size at index j
      absorption += beta(nj) * c1 * _u[_qp](j);
      emission += (j == 1 ? 2.0 : 1.0) * alpha(nj) * _u[_qp](j);
    }
    residual(0) = -(_generation - _sink * c1 - absorption + emission);
  }

  // Components i >= 1: cluster size n = i+1 >= 2
  // dC_n/dt = growth_in - growth_out + emit_in - emit_out
  // where:
  //   This form follows the single-species cluster flux definition and preserves mass:
  //   dC_n/dt = J_{n-1->n} - J_{n->n+1}, with
  //   J_{n->n+1} = beta(n)*C_1*C_n - alpha(n+1)*C_{n+1}. For n=2,
  //   J_{1->2} = beta(1)*C_1^2 - alpha(2)*C_2, so the absorption term is
  //   beta(1)*C_1^2 with no 1/2 factor.
  //   growth_in  = beta(n-1)*C_1*C_{n-1}  [cluster of size n-1 absorbs a monomer]
  //              = beta(1)*C_1^2          [special case n=2: monomer is also C_{n-1}]
  //   growth_out = beta(n)*C_1*C_n         [cluster of size n absorbs a monomer]
  //   emit_in    = alpha(n+1)*C_{n+1}      [cluster of size n+1 emits a monomer]
  //   emit_out   = alpha(n)*C_n            [cluster of size n emits a monomer]
  for (auto i = 1; i < n_comp; ++i)
  {
    const unsigned int n = i + 1;
    const auto c_n = _u[_qp](i);
    const auto c_nm1 = _u[_qp](i - 1); // for i=1: c_nm1 = c(0) = monomer

    const auto growth_in = beta(n - 1) * c1 * c_nm1;
    const auto growth_out = beta(n) * c1 * c_n;
    const auto emit_in = (i + 1 < n_comp) ? alpha(n + 1) * _u[_qp](i + 1) : 0.0;
    const auto emit_out = alpha(n) * c_n;

    residual(i) = -(growth_in - growth_out + emit_in - emit_out);
  }
}

// ---------------------------------------------------------------------------
// Non-AD Jacobian: intra-variable coupling via setJacobian(row, col, value)
// row = residual component index (cluster size n = row+1)
// col = variable component being differentiated with respect to (size n = col+1)
// ---------------------------------------------------------------------------
template <>
void
GenericClusterDynamicsNodalKernelTempl<false>::computeQpJacobian()
{
  const auto n_comp = static_cast<unsigned int>(_u[_qp].size());
  const Real c1 = _u[_qp](0);

  // --- Row 0: monomer ---
  // F(0) = -(G_1 - k_s*c(0) - 2*beta(1)*c(0)^2
  //          - sum_{j>=1} beta(j+1)*c(0)*c(j)
  //          + 2*alpha(2)*c(1)
  //          + sum_{j>=2} alpha(j+1)*c(j))

  // d F(0)/d c(0) = k_s + 4*beta(1)*c(0) + sum_{j>=1} beta(j+1)*c(j)
  {
    Real d00 = _sink + 4.0 * beta(1) * c1;
    for (unsigned int j = 1; j < n_comp; ++j)
      d00 += beta(j + 1) * _u[_qp](j);
    setJacobian(0, 0, d00);
  }

  // d F(0)/d c(j) for j >= 1: beta(j+1)*c(0) - mu(j)*alpha(j+1)
  // where mu(j) = 2 for j=1 (dimer dissociation releases 2 monomers), 1 otherwise
  for (unsigned int j = 1; j < n_comp; ++j)
  {
    const Real mu = (j == 1) ? 2.0 : 1.0;
    setJacobian(0, j, beta(j + 1) * c1 - mu * alpha(j + 1));
  }

  // --- Rows i >= 1: cluster of size n = i+1 ---
  for (unsigned int i = 1; i < n_comp; ++i)
  {
    const unsigned int n = i + 1;

    // d F(i)/d c(0): coupling to monomer
    // The n=2 case keeps a separate Jacobian entry because growth_in = beta(1)*C_1^2,
    // so differentiating with respect to the monomer concentration gives
    // 2*beta(1)*C_1 rather than beta(n-1)*C_{n-1}.
    // n=2: F(1) = -(beta(1)*c(0)^2 - beta(2)*c(0)*c(1) + ...)
    //   -> d F(1)/d c(0) = -(2*beta(1)*c(0) - beta(2)*c(1))
    // n>2: F(i) = -(beta(n-1)*c(0)*c(i-1) - beta(n)*c(0)*c(i) + ...)
    //   -> d F(i)/d c(0) = -(beta(n-1)*c(i-1) - beta(n)*c(i))
    if (n == 2)
      setJacobian(i, 0, -(2.0 * beta(1) * c1 - beta(n) * _u[_qp](i)));
    else
      setJacobian(i, 0, -(beta(n - 1) * _u[_qp](i - 1) - beta(n) * _u[_qp](i)));

    // d F(i)/d c(i-1): coupling to the next-smaller cluster (for n > 2 only;
    // for n=2 this is c(0) which is already handled above)
    if (i > 1)
      setJacobian(i, i - 1, -beta(n - 1) * c1);

    // d F(i)/d c(i): self (diagonal)
    // F(i) = -(... - beta(n)*c(0)*c(i) ... - alpha(n)*c(i))
    // -> d F(i)/d c(i) = beta(n)*c(0) + alpha(n)
    setJacobian(i, i, beta(n) * c1 + alpha(n));

    // d F(i)/d c(i+1): coupling to the next-larger cluster
    // F(i) = -(... + alpha(n+1)*c(i+1) ...)
    // -> d F(i)/d c(i+1) = -alpha(n+1)
    if (i + 1 < n_comp)
      setJacobian(i, i + 1, -alpha(n + 1));
  }
}

template <>
void
GenericClusterDynamicsNodalKernelTempl<true>::computeQpJacobian()
{
  mooseError("Internal error: computeQpJacobian should never be called for the AD version");
}

template class GenericClusterDynamicsNodalKernelTempl<false>;
template class GenericClusterDynamicsNodalKernelTempl<true>;
