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

class FirstProblemPlane;

class Snapper
{
  public:
    Snapper(Adapt* a, Tag* st, bool is);
    bool setVert(Entity* v, apf::CavityOp* o);
    bool run();
    bool dug;
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

class FirstProblemPlane
{
  public:
    FirstProblemPlane(Adapt* a, Tag* st);
    void setVertex(Entity* v);
    void setBadElements(apf::Up& badElements);
    void getCandidateEdges(std::vector<Entity*> &edges);
  private:
    Adapt* adapter;
    Tag* snapTag;
    Entity* vert;
    apf::Up problemRegions;
    Entity* problemFace;
    Entity* problemRegion;
    Vector intersection;
    apf::Up commEdges;
    double tol;
    bool find();
    void findCandidateEdges(std::vector<Entity*> &edges);
    bool intersectRayFace(const Ray& ray, const std::vector<Vector>& coords,
    	Vector& intersection, bool& isInf);
    void findCommonEdges(apf::Up& cpRegions);
};

Entity* getTetFaceOppositeVert(Mesh* m, Entity* e, Entity* v);
void getFaceCoords(Mesh* m, Entity* face, std::vector<Vector>& coords);
Vector getCenter(Mesh* mesh, Entity* face);
bool isLowInHigh(Mesh* mesh, Entity* highEnt, Entity* lowEnt);

}

#endif
