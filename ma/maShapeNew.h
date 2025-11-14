#ifndef MA_SHAPE_NEW
#define MA_SHAPE_NEW

#include "maMesh.h"
#include "maTables.h"
#include "maOperator.h"
#include "maCollapse.h"
#include "maEdgeSwap.h"
#include "maDoubleSplitCollapse.h"
#include "maSingleSplitCollapse.h"
#include "maFaceSplitCollapse.h"
#include "maFaceSwap.h"

namespace ma {
class Adapt;

class FixShape : public Operator
{
  public:
  Adapt* a;
  Mesh* mesh;
  Collapse collapse;
  SingleSplitCollapse splitCollapse;
  DoubleSplitCollapse doubleSplitCollapse;
  FaceSplitCollapse faceSplitCollapse;
  EdgeSwap* edgeSwap;
  Splits split;
  Entity* badTet;

  int numCollapse=0;
  int numEdgeSwap=0;
  int numFaceSwap=0;
  int numSplitReposition=0;
  int numRegionCollapse=0;
  int numEdgeSplitCollapse=0;
  int numFaceSplitCollapse=0;
  int numDoubleSplitCollapse=0;

  int numOneShortEdge=0;
  int numTwoShortEdge=0;
  int numThreeShortEdge=0;
  int numMoreShortEdge=0;
  int numOneLargeAngle=0;
  int numTwoLargeAngles=0;
  int numThreeLargeAngles=0;

  FixShape(Adapt* adapt);
  void apply();

  int getTargetDimension();
  bool shouldApply(Entity* e);
  bool requestLocality(apf::CavityOp* o);

  bool collapseEdge(Entity* edge);
  bool collapseToAdjacent(Entity* edge);
  double getWorstShape(EntityArray& tets, Entity*& worst);
  Entity* getLongestEdge(Entity* edges[3]);

  Vector avgCavityPos(Entity* vert);
  void repositionVertex(Entity* vert);
  bool splitReposition(Entity* edge);
  bool collapseRegion(Entity* tet, Entity* problemEnts[4]);
  
  bool isShortEdge(Entity* tet);
  bool isOneLargeAngle(Entity* tet, Entity*& worstTriangle);
  bool isTwoLargeAngles(Entity* tet, Entity* problemEnts[4]);
  
  bool fixShortEdge(Entity* tet);
  bool fixOneLargeAngle(Entity* tet);
  void fixTwoLargeAngles(Entity* tet, Entity* problemEnts[4]);
  void fixThreeLargeAngles(Entity* tet, Entity* problemEnts[4]);

  void resetCounters();
  int collect(int val);
  void printNumOperations();
  void printBadTypes();
  void printBadShape(Entity* badTet);
  void printNumTypes();
};

void fixElementShapesNew(Adapt* a);
}

#endif