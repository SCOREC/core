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
  EntityVector newEdges;
  EntityVector newTriangles;
  EntitySet elementsToChange;
  std::map<Entity*,Entity*> old2New;

  bool checkTopo(Entity* edge);
  void collectEdges(Entity* edge);
  double collectOldTets();
  void computeElementSets();
  void rebuildEdges();
  void rebuildTriangles();
  void rebuildQuads();
  void rebuildTets();
  void rebuildPrisms();
  void rebuildElements();
  bool checkValidity(double oldMinQuality);
  void destroyNewElements();
  void destroyOldElements();
  void printInfo();
};

}

#endif
