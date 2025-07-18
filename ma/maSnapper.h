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

namespace apf {
class CavityOp;
}

namespace ma {

/* tries to make room for a vertex to snap into
   a mesh by collapsing edges in the direction
   of desired snapping */

class FirstProblemPlane;

class Snapper
{
  public:
    int numFailed = 0;
    int numSnapped = 0;
    int numCollapseToVtx = 0;
    int numCollapse = 0;
    int numSwap = 0;
    int numSplitCollapse = 0;

    Snapper(Adapt* a, Tag* st, bool is);
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
    EdgeSwap* edgeSwap;

    bool tryCollapseToVertex(FirstProblemPlane* FPP);
    bool tryCollapseTetEdges(FirstProblemPlane* FPP);
    bool tryReduceCommonEdges(FirstProblemPlane* FPP);
    bool trySwapOrSplit(FirstProblemPlane* FPP);
};

struct Ray{
  Vector start;
  Vector dir;
};

class FirstProblemPlane
{
  public:
    FirstProblemPlane(Adapt* a, Tag* st);
    void setVertex(Entity* v);
    void setBadElements(apf::Up& badElements);
    void getCandidateEdges(std::vector<Entity*> &edges);
    Entity* vert;
    Entity* problemFace;
    apf::Up commEdges;
    apf::Up problemRegions;
    Entity* problemRegion;
    Tag* snapTag;
  private:
    Adapt* adapter;
    Vector intersection;
    double tol;
    bool find();
    void findCandidateEdges(std::vector<Entity*> &edges);
    bool intersectRayFace(const Ray& ray, const std::vector<Vector>& coords,
    	Vector& intersection, bool& isInf);
    void findCommonEdges(apf::Up& cpRegions);
};

int getTetStats(Adapt* a, Entity* vert, Entity* face, Entity* region, Entity* ents[4], double area[4]);
Entity* getTetFaceOppositeVert(Mesh* m, Entity* e, Entity* v);
void getFaceCoords(Mesh* m, Entity* face, std::vector<Vector>& coords);
Vector getCenter(Mesh* mesh, Entity* face);
bool isLowInHigh(Mesh* mesh, Entity* highEnt, Entity* lowEnt);

}

#endif
