#ifndef MA_LAYER_H
#define MA_LAYER_H

#include "maMesh.h"

namespace ma {

class Adapt;
class Refine;

void resetLayer(Adapt* a);

void findLayerBase(Adapt* a);

void allowSplitCollapseOutsideLayer(Adapt* a);

void allowSplitInLayer(Adapt* a);
void collectForLayerRefine(Refine* r);

int getDiagonalFromFlag(Adapt* a, Entity* e);
void tetrahedronize(Adapt* a);

void snapLayer(Adapt* a, Tag* snapTag);

}

#endif
