#ifndef MA_LAYER_H
#define MA_LAYER_H

#include "maMesh.h"

namespace ma {

class Adapt;
class Refine;

void initLayer(Adapt* a);
void preventChangesToLayer(Adapt* a);
void allowSplitCollapseOutsideLayer(Adapt* a);
void allowSplitInLayer(Adapt* a);
void turnLayerToTets(Adapt* a);

void collectForLayerRefine(Refine* r);
void flagNewLayerEntities(Refine* r);

int getDiagonalFromFlag(Adapt* a, Entity* e);

}

#endif
