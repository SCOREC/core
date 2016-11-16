/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_SNAPPER_H
#define MA_SNAPPER_H

#include "maCollapse.h"

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
    Snapper(Adapt* a, Tag* st, bool is);
    bool setVert(Entity* v, apf::CavityOp* o);
    bool run();
    bool dug;
    bool toFPP;
  private:
    Adapt* adapter;
    Tag* snapTag;
    Entity* vert;
    Collapse collapse;
    bool isSimple;
};

struct Ray{
  Vector start;
  Vector dir;
};

class FPPSnapper
{
  public:
    FPPSnapper(Adapt* a, Collapse& c, Tag* st, Entity* v, apf::Up& badElements);
    bool findFPP();
    bool snapToFPP();
  private:
    Adapt* adapter;
    Tag* snapTag;
    Entity* vert;
    Collapse collapse;
    apf::Up problemRegions;
    Entity* problemFace;
    Entity* problemRegion;
    Vector intersection;
    apf::Up commEdges;
    double tol;
    Entity* faceOppositeOfVert(Entity* e, Entity* v);
    void getFaceCoords(Entity* face, std::vector<Vector>& coords);
    bool intersectRayFace(const Ray& ray, const std::vector<Vector>& coords,
    	Vector& intersection, bool& isInf);
    void findCommonEdges(apf::Up& cpRegions);
};

Vector getCenter(Mesh* mesh, Entity* face);
bool isLowInHigh(Mesh* mesh, Entity* highEnt, Entity* lowEnt);

}

#endif
