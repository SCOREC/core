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
      orientElementVertices(mesh, Mesh::TET, tetVerts);
      newLayer.push_back(buildElement(a, mesh->toModel(oldTets[j]), Mesh::TET, tetVerts));
    }
  }

  static Entity* getPrismPartner(Mesh* m, Entity* prism, Entity* v)
  {
    Entity* pv[6];
    m->getDownward(prism, 0, pv);
    int i = apf::findIn(pv, 6, v);
    PCU_ALWAYS_ASSERT(i >= 0);
    return pv[(i + 3) % 6];
  }

  static Entity* findPrism(Mesh* m, Entity* edge, Entity* next, Entity* v)
  {
    EntityArray prisms;
    m->getAdjacent(edge, 3, prisms);
    for (size_t i = 0; i < prisms.getSize(); ++i) {
      if (m->getType(prisms[i]) != Mesh::PRISM) continue;
      Entity* pe[9];
      Entity* pv[6];
      m->getDownward(prisms[i], 1, pe);
      m->getDownward(prisms[i], 0, pv);
      if (apf::findIn(pe, 9, next) >= 0 && apf::findIn(pv, 6, v) >= 0)
        return prisms[i];
    }
    return 0;
  }

  static void buildBoundaryTris(Adapt* a, Entity* edge, Entity* va, Entity* vb, Entity* vc, Entity* vd)
  {
    Mesh* m = a->mesh;
    Model* c = m->toModel(edge);
    if (m->getModelType(c) != 2) return;
    Entity* triA[3]{vc, vd, va};
    Entity* triB[3]{vc, vd, vb};
    Entity* tris[2];
    tris[0] = buildElement(a, c, Mesh::TRIANGLE, triA);
    tris[1] = buildElement(a, c, Mesh::TRIANGLE, triB);
    if (getFlag(a, edge, LAYER_BASE))
      for (int i = 0; i < 2; ++i)
        setFlagOnClosure(a, tris[i], LAYER_BASE);
  }

  void LayerEdgeSwap::rebuildPrisms()
  {
    Entity* va = oldVerts[0];
    Entity* vb = oldVerts[1];
    Entity* vc = newVerts[0];
    Entity* vd = newVerts[1];
    buildBoundaryTris(a, edges[0], va, vb, vc, vd);
    for (size_t e = 0; e + 1 < edges.size(); ++e) {
      Entity* prismC = findPrism(mesh, edges[e], edges[e+1], vc);
      Entity* prismD = findPrism(mesh, edges[e], edges[e+1], vd);
      PCU_ALWAYS_ASSERT(prismC && prismD);
      Entity* va2 = getPrismPartner(mesh, prismC, va);
      Entity* vb2 = getPrismPartner(mesh, prismC, vb);
      Entity* vc2 = getPrismPartner(mesh, prismC, vc);
      Entity* vd2 = getPrismPartner(mesh, prismD, vd);
      buildBoundaryTris(a, edges[e+1], va2, vb2, vc2, vd2);

      Entity* vertsA[6]{vc, vd, va, vc2, vd2, va2};
      Entity* vertsB[6]{vc, vd, vb, vc2, vd2, vb2};
      orientElementVertices(mesh, Mesh::PRISM, vertsA);
      orientElementVertices(mesh, Mesh::PRISM, vertsB);
      newLayer.push_back(buildElement(a, mesh->toModel(prismC), Mesh::PRISM, vertsA));
      setFlagOnClosure(a, newLayer.back(), LAYER);
      newLayer.push_back(buildElement(a, mesh->toModel(prismD), Mesh::PRISM, vertsB));
      setFlagOnClosure(a, newLayer.back(), LAYER);

      va = va2; vb = vb2; vc = vc2; vd = vd2;
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
      if (newLayer.size() > 0) ma_dbg::flagEntity(a, 3, "new_layer", &newLayer[0], newLayer.size());

      apf::writeVtkFiles("mesh_tets", mesh, 3);
      apf::writeVtkFiles("mesh_faces", mesh, 2);
      apf::writeVtkFiles("mesh_edges", mesh, 1);
    });
    printf("COLLECTED\n");
    exit(1);
  }
}