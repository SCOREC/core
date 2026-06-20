/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_SNAPPER_H
#define MA_SNAPPER_H

#include "maCollapse.h"
#include "maSingleSplitCollapse.h"
#include "maDoubleSplitCollapse.h"
#include "maEdgeSwap.h"
#include "maReposition.h"
#include "maFirstProblemPlane.h"

namespace apf {
class CavityOp;
}

namespace ma {

/* tries to make room for a vertex to snap into
   a mesh by collapsing edges in the direction
   of desired snapping */

class Snapper
{
  public:
    int numFailed = 0;
    int numSnapped = 0;
    int numCollapseToVtx = 0;
    int numCollapse = 0;
    int numSwap = 0;
    int numSplitCollapse = 0;

    Snapper(Adapt* a, Tag* st);
    ~Snapper();
    void setVert(Entity* v);
    Entity* getVert();
    bool requestLocality(apf::CavityOp* o);
    bool trySimpleSnap();
    bool run();
  private:
    Adapt* adapt;
    Mesh* mesh;
    Tag* snapTag;
    Entity* vert;
    Collapse collapse;
    SingleSplitCollapse splitCollapse;
    DoubleSplitCollapse doubleSplitCollapse;
    RepositionVertex reposition;
    EdgeSwap* edgeSwap;

    bool tryCollapseToVertex(FirstProblemPlane* FPP);
    bool tryCollapseTetEdges(FirstProblemPlane* FPP);
    bool tryReduceCommonEdges(FirstProblemPlane* FPP);
    bool trySwapOrSplit(FirstProblemPlane* FPP);
};

EntitySet getNextLayer(Adapt* a, EntitySet& tets);

}

#endif
