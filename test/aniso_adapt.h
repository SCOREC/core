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

void refineSnapTest(ma::Mesh* m, double sizeFactor1, double sizeFactor2)
{
  m->verify();
  apf::writeVtkFiles("before_refine_snap",m);
  AnIso sf(m, sizeFactor1, sizeFactor2);
  ma::Input* in = ma::makeAdvanced(ma::configure(m, &sf));
  ma::Adapt* a = new ma::Adapt(in);
  for (int i = 0; i < in->maximumIterations; ++i)
  {
    ma::coarsen(a);
    ma::refine(a);
    ma::snap(a);
  }
  m->verify();
  apf::writeVtkFiles("after_refine_snap",m);
}