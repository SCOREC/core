#include "maAdapt.h"
#include "maLayerEdgeSwap.h"
#include "maLayer.h"
#include "maCrawler.h"
#include "pcu_util.h"
#include "maDBG.h"
#include "maFirstProblemPlane.h"

namespace ma {
  LayerEdgeSwap::LayerEdgeSwap(Adapt* adapt)
  {
    a = adapt;
    mesh = a->mesh;
  }

  bool LayerEdgeSwap::run(Entity* edge)
  {
    if (!a->hasLayer) return false;
    if (!a->input->shouldSwapLayer) return false;
    if (isOnModelEdge(mesh, edge)) return false;
    if (!getFlag(a, edge, LAYER)) return false;
    if (!checkTopo(edge)) return false;

    resetState();

    collectEdges(edge);
    computeElementSets();
    double oldMinQuality = collectOldTets();
    rebuildTets();
    rebuildPrisms();

    bool valid = checkValidity(oldMinQuality);
    if (valid) destroyOldElements();
    else destroyNewElements();
    
    return valid;
  }

  void LayerEdgeSwap::resetState()
  {
    edges.clear();
    oldTets.clear();
    oldVerts.clear();
    newVerts.clear();
    newLayer.clear();
    elementsToChange.clear();
  }

  bool LayerEdgeSwap::checkTopo(Entity* edge)
  {
    int nTets = 0;
    EntityArray adjElems;
    mesh->getAdjacent(edge, 3, adjElems);
    for (size_t i=0; i<adjElems.getSize(); ++i)
    if (mesh->getType(adjElems[i]) == Mesh::TET)
      nTets++;
    return nTets == 2;
  }

  void LayerEdgeSwap::collectEdges(Entity* edge)
  {
    Entity* verts[2];
    HasFlag p(a, CHECKED);
    while (true) {
      edges.push_back(edge);
      mesh->getDownward(edge, 0, verts);
      oldVerts.push_back(verts[0]);
      oldVerts.push_back(verts[1]);
      setFlag(a, edge, CHECKED);

      EntityArray faces;
      mesh->getAdjacent(edge, 2, faces);
      for (size_t i = 0; i < faces.getSize(); ++i)
        if (getFlag(a, faces[i], LAYER) && mesh->getType(faces[i]) == Mesh::TRIANGLE){
          newVerts.push_back(getTriVertOppositeEdge(mesh,faces[i],edge));
        }
      
      edge = getOtherEdge(mesh, edge, p);
      if (!edge) break;
    }
    for (int i=0; i<edges.size(); i++) clearFlag(a, edges[i], CHECKED);
  }

  void LayerEdgeSwap::computeElementSets()
  {
    int dim = mesh->getDimension();
    for (size_t i = 0; i < edges.size(); ++i) {
      EntityArray elements;
      mesh->getAdjacent(edges[i], dim, elements);
      for (size_t j = 0; j < elements.getSize(); ++j)
        elementsToChange.insert(elements[j]);
    }
  }

  double LayerEdgeSwap::collectOldTets()
  {
    EntityArray adjacentElements;
    mesh->getAdjacent(edges[0], 3, adjacentElements);
    double minQuality = 0.0;
    for (size_t i = 0; i < adjacentElements.getSize(); ++i)
      if (mesh->getType(adjacentElements[i]) == Mesh::TET) {
        oldTets.push_back(adjacentElements[i]);
        double quality = getAndCacheQuality(a, adjacentElements[i]);
        if (oldTets.size() == 1 || quality < minQuality)
          minQuality = quality;
      }
    PCU_ALWAYS_ASSERT(oldTets.size() == 2);
    return minQuality;
  }

  void LayerEdgeSwap::rebuildTets() 
  {
    Entity* oldFace = getTetFaceOppositeVert(mesh, oldTets[0], newVerts[0]);
    Entity* oppVert = getTriVertOppositeEdge(mesh, oldFace, edges[0]);
    for (size_t j=0; j<2; ++j) {
      Entity* tetVerts[4]{newVerts[0], newVerts[1], oppVert, oldVerts[j]};
      newLayer.push_back(buildElement(a, mesh->toModel(oldTets[j]), Mesh::TET, tetVerts));
    }
  }

  void LayerEdgeSwap::rebuildPrisms()
  {
    int nv=0;
    int ov=0;
    for (size_t e = 1; e < edges.size(); e+=2) {
      EntityArray prisms;
      mesh->getAdjacent(edges[e], 3, prisms);
      for (size_t p = 0; p < prisms.size(); p++) {
        if (mesh->getType(prisms[p]) == Mesh::TET) continue;
        Entity* verts[6]{newVerts[nv], newVerts[nv+1], newVerts[nv+2], newVerts[nv+3], oldVerts[ov], oldVerts[ov+2]};
        newLayer.push_back(buildElement(a, mesh->toModel(prisms[p]), Mesh::PRISM, verts));
        setFlagOnClosure(a, newLayer.back(), LAYER);
        ov++;
      }
      nv+=2;
    }
  }

  bool LayerEdgeSwap::checkValidity(double oldMinQuality)
  {
    for (size_t i = 0; i < newLayer.size(); ++i) {
      if (mesh->getType(newLayer[i]) == Mesh::TET) {
        if (getAndCacheQuality(a, newLayer[i]) < oldMinQuality)
          return false;
      } 
      else if (!isLayerElementOk(mesh, newLayer[i]))
        return false;
    }
    return true;
  }

  void LayerEdgeSwap::destroyNewElements()
  {
    for (size_t i = 0; i < newLayer.size(); ++i)
      destroyElement(a, newLayer[i]);
  }

  void LayerEdgeSwap::destroyOldElements()
  {
    APF_ITERATE(EntitySet, elementsToChange, it)
      destroyElement(a, *it);
  }
  
  void LayerEdgeSwap::printInfo()
  {
    ma_dbg::useFieldInfo(a, [this] {
      ma_dbg::flagEntity(a, 1, "edges_to_swap", &edges[0], edges.size());
      ma_dbg::flagEntity(a, 0, "old_verts", &oldVerts[0], oldVerts.size());
      ma_dbg::flagEntity(a, 0, "target_verts", &newVerts[0], newVerts.size());
      ma_dbg::flagEntity(a, 3, "tets_to_swap", elementsToChange);

      apf::writeVtkFiles("mesh_tets", mesh, 3);
      apf::writeVtkFiles("mesh_faces", mesh, 2);
      apf::writeVtkFiles("mesh_edges", mesh, 1);
    });
    printf("COLLECTED\n");
    exit(1);
  }
}