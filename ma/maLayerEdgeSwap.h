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
  EntityVector oldVerts;
  EntityVector newVerts;
  EntityVector newLayer;
  EntityVector newEdges;
  EntityVector newTriangles;
  std::map<Entity*,Entity*> old2New;
  EntitySet elementsToChange;

  void collectEdges(Entity* edge);
  void computeElementSets();
  void rebuildEdges();
  void rebuildTriangles();
  void rebuildQuads();
  void rebuildTets();
  void rebuildElements();
  void destroyOldElements();
  void printInfo();
};

}

#endif
