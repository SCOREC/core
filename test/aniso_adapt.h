#ifndef TEST_ANISO_ADAPT_H
#define TEST_ANISO_ADAPT_H

#include <iostream>
#include <lionPrint.h>
#include <pcu_util.h>
#include <gmi_mesh.h>
#include <apf.h>
#include <apfMDS.h>
#include <ma.h>
#include <maCoarsen.h>
#include <maAdapt.h>
#include <maRefine.h>
#include <maShape.h>
#include <maSnap.h>
#include <apfGeometry.h>
#include <functional>
/*
 Test some of the individual components in mesh adaptation to make sure that they are 
 functioning as intended. Right now it only tests coarsen refinement and snapping but
 can be expanded to test more in the future. This has been tested with mds, simmetrix,
 and capstone meshes.
*/
class AnIso : public ma::AnisotropicFunction
{
  public:
    AnIso(ma::Mesh* m, double sf1, double sf2) :
      sizeFactor1(sf1), sizeFactor2(sf2)
    {
      average = ma::getAverageEdgeLength(m);
      ma::getBoundingBox(m, lower, upper);
    }
    virtual void getValue(ma::Entity*, ma::Matrix& R, ma::Vector& H)
    {
      double h = average/sizeFactor1;
      H = ma::Vector(h, h, h/sizeFactor2);
      R = ma::Matrix(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0
      );
    }
  private:
    double sizeFactor1, sizeFactor2, average;
    ma::Vector lower, upper;
};

void measureQuality(ma::Mesh* m, double& avgQuality, double& minQuality)
{
  int d = m->getDimension();
  apf::MeshEntity* elem;
  apf::MeshIterator* it = m->begin(d);
  ma::IdentitySizeField I(m);
  while ((elem = m->iterate(it))) {
    double q = ma::measureElementQuality(m, &I, elem);
    avgQuality += q;
    minQuality = std::min(minQuality, q);
  }
  m->end(it);
  avgQuality = avgQuality / m->count(d);
}

//make sure refinement is done
bool noLongEdges(ma::Adapt* a)
{
  apf::MeshEntity* edge;
  apf::MeshIterator* it = a->mesh->begin(1);
  while ((edge = a->mesh->iterate(it)))
    if (a->sizeField->measure(edge) > ma::MAXLENGTH)
      return false;
  a->mesh->end(it);
  return true;
}

int countEdges(ma::Mesh* m)
{
  return m->count(1);
}

ma::Mesh* coarsenForced(ma::Mesh* m)
{
  m->verify();
  AnIso sf(m, .5, 1);
  ma::Input* in = ma::makeAdvanced(ma::configure(m, &sf));
  in->shouldForceAdaptation = true;
  ma::Adapt* a = new ma::Adapt(in);
  double avgQualBefore, avgQualAfter, minQualBefore, minQualAfter;

  measureQuality(m, avgQualBefore, minQualBefore);
  double averageBefore = ma::getAverageEdgeLength(m);
  int edgesBefore = countEdges(m);

  // ma::coarsenMultiple(a);

  measureQuality(m, avgQualAfter, minQualAfter);

  //make sure that coarsening is happening and it doesn't make the mesh invalid
  PCU_ALWAYS_ASSERT(edgesBefore > countEdges(m));
  PCU_ALWAYS_ASSERT(averageBefore < ma::getAverageEdgeLength(m));
  PCU_ALWAYS_ASSERT(minQualAfter >= in->validQuality);

  m->verify();
  delete a;
  // cleanup input object and associated sizefield and solutiontransfer objects
  if (in->ownsSizeField)
    delete in->sizeField;
  if (in->ownsSolutionTransfer)
    delete in->solutionTransfer;
  delete in;
  return m;
}

ma::Mesh* coarsenRegular(ma::Mesh* m)
{
  m->verify();
  AnIso sf(m, .5, 1);
  ma::Input* in = ma::makeAdvanced(ma::configure(m, &sf));
  ma::Adapt* a = new ma::Adapt(in);
  double avgQualBefore, avgQualAfter, minQualBefore, minQualAfter;

  measureQuality(m, avgQualBefore, minQualBefore);
  double averageBefore = ma::getAverageEdgeLength(m);
  int edgesBefore = countEdges(m);

  // ma::coarsenMultiple(a);

  measureQuality(m, avgQualAfter, minQualAfter);

  //Makes sure that coarsening is happening and it isn't decreasing the quality of the mesh
  PCU_ALWAYS_ASSERT(edgesBefore > countEdges(m));
  PCU_ALWAYS_ASSERT(averageBefore < ma::getAverageEdgeLength(m));
  PCU_ALWAYS_ASSERT(fabs(minQualBefore - minQualAfter) < 0.001f || minQualBefore < minQualAfter);

  m->verify();
  delete a;
  // cleanup input object and associated sizefield and solutiontransfer objects
  if (in->ownsSizeField)
    delete in->sizeField;
  if (in->ownsSolutionTransfer)
    delete in->solutionTransfer;
  delete in;
  return m;
}

//Makes sure that all of the vertices snapped to the model
bool allVertsOnModel(ma::Adapt* a)
{
  if (!a->input->shouldSnap)
    return true;
  ma::Mesh* m = a->mesh;
  int dim = m->getDimension();
  ma::Iterator* it = m->begin(0);
  ma::Entity* vert;
  while ((vert = m->iterate(it))) {
    int md = m->getModelType(m->toModel(vert));
    if (dim == 3 && md == 3)
      continue;
    ma::Vector modelPoint, param;
    m->getPoint(vert,0,modelPoint);
    m->getParam(vert,param);
    ma::Model* model = m->toModel(vert);
    m->snapToModel(model,param,modelPoint);
    ma::Vector position = ma::getPosition(m, vert);
    if (apf::areClose(modelPoint, position, 1e-12))
      continue;
  }
  m->end(it);
  return true;
}

ma::Mesh* refineSnapTest(ma::Mesh* m)
{
  m->verify();
  AnIso sf(m, 3, 1);
  ma::Input* in = ma::makeAdvanced(ma::configure(m, &sf));
  ma::Adapt* a = new ma::Adapt(in);
  int edgesBefore = countEdges(m);
  double averageBefore = ma::getAverageEdgeLength(m);

  for (int i = 0; i < in->maximumIterations; ++i)
  {
    ma::refine(a);
    ma::snap(a);
  }

  double avgQualAfter, minQualAfter;
  measureQuality(m, avgQualAfter, minQualAfter);

  //Makes sure that refinement is happening
  PCU_ALWAYS_ASSERT(edgesBefore < countEdges(m));
  PCU_ALWAYS_ASSERT(averageBefore > ma::getAverageEdgeLength(m));

  PCU_ALWAYS_ASSERT(noLongEdges(a));
  PCU_ALWAYS_ASSERT(allVertsOnModel(a));
  PCU_ALWAYS_ASSERT(minQualAfter >= 0);

  m->verify();
  delete a;
  // cleanup input object and associated sizefield and solutiontransfer objects
  if (in->ownsSizeField)
    delete in->sizeField;
  if (in->ownsSolutionTransfer)
    delete in->solutionTransfer;
  delete in;
  return m;
}

ma::Mesh* shapeTest(ma::Mesh* m)
{
  m->verify();
  AnIso sf(m, 8, 1);
  ma::Input* in = ma::makeAdvanced(ma::configure(m, &sf));
  adapt(in);
  delete in;
  return m;
}

void adaptTests(const std::function<ma::Mesh*()>& createMesh)
{
  //Mesh created multiple times to compare adaptation
  ma::Mesh* meshReg = createMesh();
  apf::writeVtkFiles("startMesh", meshReg);

  shapeTest(meshReg);
  apf::writeVtkFiles("afterShape", meshReg);

  // refineSnapTest(meshReg);
  // apf::writeVtkFiles("afterRefine", meshReg);

  // coarsenRegular(meshReg);
  // apf::writeVtkFiles("afterCoarsen", meshReg);

  // ma::Mesh* meshForce = coarsenForced(refineSnapTest(createMesh()));
  // apf::writeVtkFiles("afterForcedCoarsen", meshForce);

  // //Make sure setting to force coarsen is functioning
  // PCU_ALWAYS_ASSERT(countEdges(meshReg) > countEdges(meshForce));
  // PCU_ALWAYS_ASSERT(ma::getAverageEdgeLength(meshReg) < ma::getAverageEdgeLength(meshForce));

  meshReg->destroyNative();
  apf::destroyMesh(meshReg);
}

#endif