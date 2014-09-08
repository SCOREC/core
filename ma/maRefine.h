/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_REFINE_H
#define MA_REFINE_H

#include "maMesh.h"
#include "maTables.h"

namespace ma {

class Adapt;

class Refine
{
  public:
    Refine(Adapt* a);
    ~Refine();
    Adapt* adapt;
    Tag* numberTag;
    Tag* vertPlaceTag;
    EntityArray toSplit[4];
    apf::DynamicArray<EntityArray> newEntities[4];
    bool shouldCollect[4];
};

void addEdgePreAllocation(Refine* r, Entity* edge, int counts[4]);
void allocateRefine(Refine* r, int counts[4]);
void addEdgePostAllocation(Refine* r, Entity* edge, int counts[4]);

void resetCollection(Refine* r);
void collectForTransfer(Refine* r);
void collectForMatching(Refine* r);

void transferElements(Refine* r);
void forgetNewEntities(Refine* r);
void destroySplitElements(Refine* r);
void cleanSplitVerts(Refine* r);

void splitElements(Refine* r);
void processNewElements(Refine* r);
void cleanupAfter(Refine* r);

bool refine(Adapt* a);

Entity* buildSplitElement(
    Refine* r,
    Entity* parent,
    int type,
    Entity** verts);

Entity* findSplitVert(Refine* r, int dimension, int id);
Entity* findSplitVert(Refine* r, Entity* parent);
Entity* findSplitVert(Refine* r, Entity* v0, Entity* v1);
Entity* findPlacedSplitVert(Refine* r, Entity* v0, Entity* v1, double& place);

typedef void (*SplitFunction)(Refine* r, Entity* p, Entity** v);

int matchEntityToTemplate(Adapt* a, Entity* e, Entity** vo);
int matchToTemplate(int type, Entity** vi, int code, Entity** vo);

}

#endif
