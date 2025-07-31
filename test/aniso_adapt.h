#include <iostream>
#include <cstdlib>
#include <filesystem>

#include <lionPrint.h>
#include <pcu_util.h>
#include <gmi_mesh.h>
#include <apf.h>
#include <apfMDS.h>
#include <ma.h>
#include "maCoarsen.h"
#include "maAdapt.h"
#include "maRefine.h"
#include "maShape.h"
#include "maSnap.h"

class AnIso : public ma::AnisotropicFunction
{
  public:
    AnIso(ma::Mesh* m, double sf1, double sf2) :
      mesh(m), sizeFactor1(sf1), sizeFactor2(sf2)
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
    ma::Mesh* mesh;
    double sizeFactor1, sizeFactor2, average;
    ma::Vector lower, upper;
};

void measureQuality(ma::Mesh* m, double& avgQuality, double& minQuality)
{
  int d = m->getDimension();
  apf::MeshEntity* elem;
  apf::MeshIterator* elems = m->begin(d);
  ma::IdentitySizeField I(m);
  while ((elem = m->iterate(elems))) {
    double q = ma::measureElementQuality(m, &I, elem);
    avgQuality += q;
    minQuality = fmin(minQuality, q);
  }
  m->end(elems);
  avgQuality = avgQuality / m->count(d);
}

int countEdges(ma::Mesh* m)
{
  return m->count(1);
}

void coarsenTest(ma::Mesh* m)
{
  m->verify();
  apf::writeVtkFiles("before_coarsen_test",m);
  AnIso sf(m, .5, 1);
  ma::Input* in = ma::makeAdvanced(ma::configure(m, &sf));
  ma::Adapt* a = new ma::Adapt(in);
  double avgQualBefore, avgQualAfter, minQualBefore, minQualAfter;

  measureQuality(m, avgQualBefore, minQualBefore);
  double averageBefore = ma::getAverageEdgeLength(m);
  int edgesBefore = countEdges(m);
  ma::coarsenMultiple(a);
  measureQuality(m, avgQualAfter, minQualAfter);
  double averageAfter = ma::getAverageEdgeLength(m);
  int edgesAfter = countEdges(m);

  PCU_ALWAYS_ASSERT(edgesBefore > edgesAfter);
  PCU_ALWAYS_ASSERT(averageBefore < averageAfter);
  PCU_ALWAYS_ASSERT(minQualBefore <= minQualAfter);
  // PCU_ALWAYS_ASSERT(avgQualBefore <= avgQualAfter);


  //TODO: Check that force coarsen again doesn't result in any changes
  m->verify();
  apf::printStats(m);
  apf::writeVtkFiles("after_coarsen_test",m);
  delete a;
  // cleanup input object and associated sizefield and solutiontransfer objects
  if (in->ownsSizeField)
    delete in->sizeField;
  if (in->ownsSolutionTransfer)
    delete in->solutionTransfer;
  delete in;
}

void refineSnapTest(ma::Mesh* m)
{
  m->verify();
  apf::writeVtkFiles("before_adapt",m);
  AnIso sf(m, 2, 1);
  ma::Input* in = ma::makeAdvanced(ma::configure(m, &sf));
  ma::Adapt* a = new ma::Adapt(in);
  int edgesBefore = countEdges(m);
  double averageBefore = ma::getAverageEdgeLength(m);

  for (int i = 0; i < in->maximumIterations; ++i)
  {
    ma::refine(a);
    ma::snap(a);
  }
  int edgesAfter = countEdges(m);
  double averageAfter = ma::getAverageEdgeLength(m);
  PCU_ALWAYS_ASSERT(edgesBefore < edgesAfter);
  PCU_ALWAYS_ASSERT(averageBefore > averageAfter);


  m->verify();
  apf::writeVtkFiles("after_adapt",m);
}