
#ifndef MA_LAYER_COLLAPSE_H
#define MA_LAYER_COLLAPSE_H

#include "maMesh.h"
#include "maCollapse.h"
#include <vector>

namespace ma {

typedef std::vector<Entity*> EntityVector;

class Adapt;

struct LayerCollapse
{
  LayerCollapse(Adapt* a_);
  /* return true iff the stacks were local enough to collapse */
  bool setup(Entity* edge);
  /* return true iff the collapse was successful */
  bool apply(double qualityToBeat);
  Adapt* a;
  Mesh* m;
  Collapse collapse;
  EntityVector edges;
  EntityVector curves[2];
  EntitySet elementsToCollapse;
  EntitySet elementsToKeep;
  EntityArray newSimplices;
  EntityVector newLayer;
private:
  void computeElementSets();
  bool involvesPyramids();
  bool checkIndividualCollapses();
  void rebuildElements();
  bool checkValidity(double qualityToBeat);
  void destroyOldElements();
  bool apply_(double qualityToBeat);
  bool setup_(Entity* edge);
  void unmark();
  void cancel();
};

}

#endif
