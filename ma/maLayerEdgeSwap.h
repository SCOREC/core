#ifndef MA_LAYEREDGESWAP_H
#define MA_LAYEREDGESWAP_H

#include "maMesh.h"
#include "maEdgeSwap.h"
#include <vector>

namespace ma {

typedef std::vector<Entity*> EntityVector;

class LayerEdgeSwap
{
  public:
  LayerEdgeSwap(Adapt* adapt);
  bool run(Entity* edge);
  void testAll();

  private:
  Adapt* a;
  Mesh* mesh;
  EntityVector edges;
  EntityVector oldTets;
  EntityVector oldVerts;
  EntityVector newVerts;
  EntityVector newLayer;
  EntitySet elementsToChange;

  void resetState();
  bool checkTopo(Entity* edge);
  void collectEdges(Entity* edge);
  double collectOldTets();
  void computeElementSets();
  void rebuildTets();
  void rebuildPrisms();
  bool checkValidity(double oldMinQuality);
  void destroyNewElements();
  void destroyOldElements();
  void printInfo();
};

}

#endif
