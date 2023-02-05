/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include "apfElement.h"
#include <apfCavityOp.h>
#include "ree.h"

namespace ree {


struct Sizefield {
  apf::Mesh* mesh;
  /* the polynomial order of the solution electric Nedelec field */
  int order;
  /* the input solution electric field obtained after the FEM solve */
  apf::Field* ef;
  /* the per element error field obtained after executing the implicit
     residual error estimator */
  apf::Field* errorField;
  double size_factor;
  /* target error = sum of all element erros / total number of elements */
  double target_error;
  double alpha;
  double beta;
  /* a temporary field storing desired sizes at elements */
  apf::Field* element_size;
  /* the resulting size field, recovered from the element_size field
     (using a local average recovery method much weaker than SPR) */
  apf::Field* size;
};

static void setupSizefield(
    Sizefield* sz,
    apf::Field* ef,
    apf::Field* error_field,
    int n,
    double alpha,
    double beta)
{
  sz->mesh = ef->getMesh();
  sz->order = ef->getShape()->getOrder();
  sz->errorField = error_field;
  sz->size_factor = 0;

  apf::MeshEntity* entity;
  int d = sz->mesh->getDimension();
  apf::MeshIterator* it = sz->mesh->begin(d);
  double total_error = 0.;
  while ((entity = sz->mesh->iterate(it))) {
    total_error += apf::getScalar(sz->errorField, entity, 0);
  }
  sz->mesh->end(it);
  sz->target_error = total_error / (sz->mesh->count(d)* n);

  sz->alpha = alpha;
  sz->beta = beta;
  sz->element_size = 0;
  sz->size = 0;
}

static double getCurrentSize(apf::Mesh* m, apf::MeshEntity* e)
{
  /* right now minimum edge length is the formula... */
  apf::Downward edges;
  int ne = m->getDownward(e,1,edges);
  double h = std::numeric_limits<double>::max();
  for (int i=0; i < ne; ++i)
    h = std::min(h, measure(m, edges[i]));
  return h;
}

static double getDesiredSize(Sizefield* sz, apf::MeshEntity* entity)
{
  double element_error = apf::getScalar(sz->errorField, entity, 0);
  double h = getCurrentSize(sz->mesh, entity);

  double p = sz->order;
  int d = sz->mesh->getDimension();
  sz->size_factor = pow(sz->target_error/element_error, (2. / (2.*p + d)));

  double h_new = h * sz->size_factor;
  if (h_new < sz->alpha*h) h_new = sz->alpha*h;
  if (h_new > sz->beta*h) h_new = sz->beta*h;
  return h_new;
}

static void getElementSizeField(Sizefield* sz)
{
  apf::Field* eSize = apf::createStepField(
      sz->mesh, "esize", apf::SCALAR);
  int d = sz->mesh->getDimension();
  apf::MeshEntity* entity;
  apf::MeshIterator* elements = sz->mesh->begin(d);
  while ((entity = sz->mesh->iterate(elements))) {
    double h = getDesiredSize(sz, entity);
    apf::setScalar(eSize, entity, 0, h);
  }
  sz->mesh->end(elements);
  sz->element_size = eSize;
}


void averageToVertex(apf::Field* ef, apf::Field* vf, apf::MeshEntity* ent)
{
  apf::Mesh* m = apf::getMesh(ef);
  apf::Adjacent elements;
  m->getAdjacent(ent, m->getDimension(), elements);
  double s=0;
  for (std::size_t i=0; i < elements.getSize(); ++i)
    s += apf::getScalar(ef, elements[i], 0);
  s /= elements.getSize();
  apf::setScalar(vf, ent, 0, s);
}

class AverageOp : public apf::CavityOp
{
public:
  AverageOp(Sizefield* sz):
    apf::CavityOp(sz->mesh),
    sizefield(sz),
    entity(0)
  {
  }
  virtual Outcome setEntity(apf::MeshEntity* e)
  {
    entity = e;
    if (apf::hasEntity(sizefield->size, entity))
      return SKIP;
    if ( !requestLocality(&entity,1))
      return REQUEST;
    return OK;
  }
  virtual void apply()
  {
    averageToVertex(sizefield->element_size,
        sizefield->size, entity);
  }
  Sizefield* sizefield;
  apf::MeshEntity* entity;
};

void averageSizeField(Sizefield* sz)
{
  sz->size = apf::createLagrangeField(sz->mesh, "size", apf::SCALAR, 1);
  AverageOp op(sz);
  op.applyToDimension(0);
}

apf::Field* getTargetEMSizeField(
    apf::Field* ef,
    apf::Field* error_field,
    int n,
    double alpha /*= 0.25*/,
    double beta /*= 2.0*/)
{
  double t0 = PCU_Time();
  Sizefield sz;
  setupSizefield(&sz, ef, error_field, n, alpha, beta);
  getElementSizeField(&sz);
  averageSizeField(&sz);
  apf::destroyField(sz.element_size);
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
    lion_eprint(1,"EM: SizeField computed in %f seconds\n",t1-t0);
  return sz.size;
}

}
