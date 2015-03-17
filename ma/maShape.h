/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_SHAPE
#define MA_SHAPE

#include "maMesh.h"
#include "maTables.h"

namespace ma {

class SizeField;
class Adapt;

double measureTriQuality(Mesh* m, SizeField* f, Entity* tri);
double measureTetQuality(Mesh* m, SizeField* f, Entity* tet);
double measureElementQuality(Mesh* m, SizeField* f, Entity* e);

double measureQuadraticTetQuality(Mesh* m, Entity* tet);

double getWorstQuality(Adapt* a, EntityArray& e);
double getWorstQuality(Adapt* a, Entity** e, size_t n);

bool isLayerElementOk(Mesh* m, Entity* e);

CodeMatch matchSliver(
    Mesh* m,
    Entity* tet);

void fixElementShapes(Adapt* a);
void printQuality(Adapt* a);

}

#endif
