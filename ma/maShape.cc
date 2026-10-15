
/******************************************************************************

  Copyright 2013 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

*******************************************************************************/
#include "maShape.h"
#include "maSize.h"
#include "maAdapt.h"
#include "maSnap.h"
#include "maOperator.h"
#include "maEdgeSwap.h"
#include "maDoubleSplitCollapse.h"
#include "maSingleSplitCollapse.h"
#include "maFaceSplitCollapse.h"
#include "maShortEdgeRemover.h"
#include "maShapeHandler.h"
#include "maBalance.h"
#include "maDBG.h"
#include <pcu_util.h>

namespace ma {

double getMinQuality(Adapt* a)
{
  PCU_ALWAYS_ASSERT(a);
  Mesh* m;
  m = a->mesh;
  PCU_ALWAYS_ASSERT(m);
  Iterator* it = m->begin(m->getDimension());
  Entity* e;
  double minqual = 1;
  while ((e = m->iterate(it))) {
    if (!apf::isSimplex(m->getType(e)))
      continue;
    double qual = a->shape->getQuality(e);
    if (qual < minqual)
      minqual = qual;
  }
  m->end(it);
  return m->getPCU()->Min<double>(minqual);
}

class QualityImprover2D : public Operator
{
  public:
    QualityImprover2D(Adapt* a)
    {
      adapter = a;
      mesh = a->mesh;
      edgeSwap = makeEdgeSwap(a);
      ns = nf = 0;
      edge = 0;
    }
    virtual ~QualityImprover2D()
    {
      delete edgeSwap;
    }
    virtual int getTargetDimension() {return 1;}
    virtual bool shouldApply(Entity* e)
    {
      if ( getFlag(adapter,e,DONT_SWAP))
        return false;
      if ( mesh->isShared(e) )
      	return false;
      edge = e;
      return true;
    }
    virtual bool requestLocality(apf::CavityOp* o)
    {
      return o->requestLocality(&edge,1);
    }
    virtual void apply()
    {
      if (edgeSwap->run(edge))
      {
	++ns;
	return;
      }
      ++nf;
    }
  private:
    Adapt* adapter;
    Mesh* mesh;
    Entity* edge;
    EdgeSwap* edgeSwap;
    int ns;
    int nf;
};

static void improveQualities2D(Adapt* a)
{
  QualityImprover2D improver(a);
  applyOperator(a, &improver);
}


double improveQualities(Adapt* a)
{
  double t0 = pcu::Time();
  if (a->mesh->getDimension() == 3)
    return 0; // TODO: implement this for 3D
  else
    improveQualities2D(a);
  double t1 = pcu::Time();
  return t1 - t0;
}


void printQuality(Adapt* a)
{
  if ( ! a->input->shouldPrintQuality)
    return;
  double minqual = cbrt(getMinQuality(a));
  print(a->mesh->getPCU(), "worst element quality is %e",  minqual);
}

}