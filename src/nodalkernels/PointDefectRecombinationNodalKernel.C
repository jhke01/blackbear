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

#include "PointDefectRecombinationNodalKernel.h"

registerMooseObject("BlackBearApp", PointDefectRecombinationNodalKernel);

InputParameters
PointDefectRecombinationNodalKernel::validParams()
{
  InputParameters params = ArrayNodalKernel::validParams();
  params.addClassDescription(
      "Point-defect recombination loss on this variable's monomer component (component 0), "
      "coupled to the other defect family's own array variable. Apply one instance per family "
      "(each coupled to the other) to represent recombination symmetrically.");
  params.addRequiredParam<Real>(
      "recombination_rate",
      "Recombination rate coefficient [m^3/s], 4*pi*r_rec*(D_i+D_v)/V_at.");
  params.addRequiredCoupledVar(
      "coupled_var", "The other defect family's grouped array variable.");
  return params;
}

PointDefectRecombinationNodalKernel::PointDefectRecombinationNodalKernel(
    const InputParameters & parameters)
  : ArrayNodalKernel(parameters),
    _recombination_rate(getParam<Real>("recombination_rate")),
    _coupled_monomer(coupledArrayValue("coupled_var")),
    _coupled_var_num(coupled("coupled_var"))
{
  if (_recombination_rate < 0.0)
    mooseError("PointDefectRecombinationNodalKernel requires recombination_rate >= 0.");
}

void
PointDefectRecombinationNodalKernel::computeQpResidual(RealEigenVector & residual)
{
  residual.resize(_count);
  residual.setZero();
  residual(0) = _recombination_rate * _u[_qp](0) * _coupled_monomer[_qp](0);
}

void
PointDefectRecombinationNodalKernel::computeQpJacobian()
{
  setJacobian(0, 0, _recombination_rate * _coupled_monomer[_qp](0));
}

void
PointDefectRecombinationNodalKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar != _coupled_var_num)
    return;

  setJacobian(0, 0, _recombination_rate * _u[_qp](0));
}
