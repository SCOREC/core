#include "maFaceSplitCollapse.h"
#include "maAdapt.h"
#include "maShape.h"
#include <pcu_util.h>

namespace ma {

  FaceSplitCollapse::FaceSplitCollapse(Adapt* a):
    faceSplit(a)
  {
    collapse.Init(a);
    oldQuality = 2;
  }

  void FaceSplitCollapse::getNewElements(EntityArray& e)
  {
    Mesh* m = getAdapt()->mesh;
    EntityArray& c = collapse.newElements;
    EntitySet& b = collapse.elementsToCollapse;
    EntityArray ab;
    m->getAdjacent(collapse.vertToKeep,3,ab);
    e.setSize(ab.getSize() - b.size() + c.getSize());
    size_t i=0;
    for (size_t j=0; j < ab.getSize(); ++j)
      if ( ! b.count(ab[j]))
	e[i++] = ab[j];
    for (size_t j=0; j < c.getSize(); ++j)
      e[i++] = c[j];
    PCU_ALWAYS_ASSERT(i==e.getSize());
  }

  bool FaceSplitCollapse::tryBothCollapses(Entity* e)
  {
    if ( ! collapse.setEdge(e))
      return false;
    if ( ! collapse.checkClass())
      return false;
    if ( ! collapse.checkTopo())
      return false;
    EntityArray& preSplit = faceSplit.getTets();
    for (size_t i=0; i < preSplit.getSize(); ++i)
      collapse.elementsToIgnore.insert(preSplit[i]);
    return collapse.tryBothDirections(oldQuality);
    collapse.elementsToIgnore.clear();
  }

  void FaceSplitCollapse::accept()
  {
    faceSplit.destroyOldElements();
    collapse.destroyOldElements();
  }

  bool FaceSplitCollapse::run(Entity* face, Entity* tet)
  {
    Adapt* a = getAdapt();
    Mesh* m = a->mesh;
    if ( ! faceSplit.setFace(face))
      return false;
    // TODO: Assert that tet & face are adjacent?
    oldQuality = getWorstQuality(a, faceSplit.getTets());
    faceSplit.makeNewElements();
    faceSplit.transfer();
    Entity* splitVerts[2];
    splitVerts[0] = faceSplit.getSplitVert();
    splitVerts[1] = getTetVertOppositeTri(m, tet, face);
    Entity* edge = findUpward(m, apf::Mesh::EDGE, splitVerts);
    if (tryBothCollapses(edge))
    {
      accept();
      return true;
    }
    else
    {
      faceSplit.cancel();
      return false;
    }
  }

  Adapt* FaceSplitCollapse::getAdapt()
  {
    return collapse.adapt;
  }

void FaceSplitCollapse::IgnoringCollapse::computeElementSets()
{
  Collapse::computeElementSets();
  APF_ITERATE(EntitySet, elementsToIgnore, it) {
    if (elementsToKeep.count(*it))
      elementsToKeep.erase(*it);
    if (elementsToCollapse.count(*it))
      elementsToCollapse.erase(*it);
  }
}

}
