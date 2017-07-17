/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include <PCU.h>
#include "crvAdapt.h"
#include "maOperator.h"
#include "maEdgeSwap.h"
#include "maDoubleSplitCollapse.h"
#include "maShortEdgeRemover.h"
#include <pcu_util.h>

#ifndef CRV_CRVSHAPEFIX_H
#define CRV_CRVSHAPEFIX_H

namespace crv {

class CrvTetFixerBase
{
  public:
    virtual void setTet(apf::MeshEntity** e) = 0;
    virtual bool requestLocality(apf::CavityOp* o) = 0;
    virtual bool run() = 0;
    virtual int getSuccessCount() = 0;
};

class CrvFaceVertFixer : public CrvTetFixerBase
{
  public:
    CrvFaceVertFixer(Adapt* a);
    ~CrvFaceVertFixer();
    virtual void setTet(apf::MeshEntity** v);
    virtual bool requestLocality(apf::CavityOp* o);
    virtual bool run();
    virtual int getSuccessCount();
  private:
    apf::Mesh2* mesh;
    apf::MeshEntity* edges[3];
    ma::EdgeSwap* edgeSwap;
    int nes;
    int nf;
};

class CrvEdgeEdgeFixer : public CrvTetFixerBase
{
  public:
    CrvEdgeEdgeFixer(Adapt* a);
    ~CrvEdgeEdgeFixer();
    virtual void setTet(apf::MeshEntity** v);
    virtual bool requestLocality(apf::CavityOp* o);
    virtual bool run();
    virtual int getSuccessCount();
  private:
    Adapt* adapter;
    apf::Mesh2* mesh;
    apf::MeshEntity* edges[2];
    ma::EdgeSwap* edgeSwap;
    ma::DoubleSplitCollapse doubleSplitCollapse;
    int nes;
    int ndsc;
    int nf;
    ma::SizeField* sf;
};

class CrvLargeAngleTetFixer : public ma::Operator
{
  public:
    CrvLargeAngleTetFixer(Adapt* a);
    virtual ~CrvLargeAngleTetFixer();
    virtual int getTargetDimension();
    enum{EDGE_EDGE, FACE_VERT};
    virtual bool shouldApply(apf::MeshEntity* e);
    virtual bool requestLocality(apf::CavityOp* o);
    virtual void apply();
    int getSuccessCount();
  private:
    Adapt* adapter;
    apf::Mesh2* mesh;
    apf::MeshEntity* tet;
    CrvEdgeEdgeFixer edgeEdgeFixer;
    CrvFaceVertFixer faceVertFixer;
    CrvTetFixerBase* fixer;
};

class CrvShortEdgeFixer : public ma::Operator
{
  public:
    CrvShortEdgeFixer(Adapt* a);
    virtual ~CrvShortEdgeFixer();
    virtual int getTargetDimension();
    virtual bool shouldApply(apf::MeshEntity* e);
    virtual bool requestLocality(apf::CavityOp* o);
    virtual void apply();
    int nr;
  private:
    Adapt* adapter;
    apf::Mesh2* mesh;
    apf::MeshEntity* element;
    ma::SizeField* sizeField;
    ma::ShortEdgeRemover remover;
    double shortEdgeRatio;
    int nf;
};

}
#endif
