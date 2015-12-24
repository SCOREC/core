/******************************************************************************

  Copyright 2015 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#include "size.h"
#include <apf.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <cassert>

namespace size {

namespace mech {

class MechanicsError : public apf::Integrator
{
  public:
    MechanicsError(Input const& in);
    void inElement(apf::MeshElement* me);
    void outElement();
    void atPoint(apf::Vector3 const& p, double dv, double w);
    apf::Field* getResidual() {return residual;}
    double getGlobalError() {return globalError;}
  private:
    apf::Mesh* mesh;
    int ndims;
    apf::Field* stress;
    apf::Field* adjoint;
    apf::Field* residual;
    apf::MeshEntity* elem;
    apf::Element* adjointElement;
    double localError;
    double globalError;
};

MechanicsError::MechanicsError(Input const& in) :
  apf::Integrator(apf::getShape(in.stress)->getOrder())
{
  mesh = apf::getMesh(in.stress);
  ndims = mesh->getDimension();
  stress = in.stress;
  adjoint = in.adjoint;
  residual = apf::createStepField(mesh, "dwr_mech_residual", apf::SCALAR);
  elem = 0;
  adjointElement = 0;
  localError = 0.0;
  globalError = 0.0;
}

void MechanicsError::inElement(apf::MeshElement* me)
{
  elem = apf::getMeshEntity(me);
  adjointElement = apf::createElement(adjoint, me);
  localError = 0.0;
}

void MechanicsError::outElement()
{
  apf::destroyElement(adjointElement);
}

void MechanicsError::atPoint(apf::Vector3 const& p, double w, double dv)
{
  apf::Matrix3x3 gradz;
  apf::getVectorGrad(adjointElement, p, gradz);
  apf::Matrix3x3 sigma;
  apf::getMatrix(stress, elem, this->ipnode, sigma);
  for (int i=0; i < ndims; ++i)
  for (int j=0; j < ndims; ++j)
    localError += sigma[i][j] * gradz[i][j] * w * dv;
  apf::setScalar(residual, elem, this->ipnode, localError);
  globalError += localError;
}

static void validateInput(Input const& in)
{
  assert(in.stress);
  assert(in.adjoint);
  assert(apf::getValueType(in.stress) == apf::MATRIX);
  assert(apf::getValueType(in.adjoint) == apf::VECTOR);
}

apf::Field* estimateError(Input const& in)
{
  validateInput(in);
  MechanicsError e(in);
  e.process(apf::getMesh(in.adjoint));
  return e.getResidual();
}

}

}
