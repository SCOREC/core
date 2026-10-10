#include "maAdapt.h"
#include "maLayerEdgeSwap.h"
#include "maLayer.h"
#include "maCrawler.h"
#include "pcu_util.h"
#include "maDBG.h"
#include "maFirstProblemPlane.h"

namespace ma {
  namespace {
    Entity* getUpward(Mesh* mesh, int type, std::initializer_list<Entity*> down)
    {
      std::vector<Entity*> entities(down.begin(), down.end());
      Entity* upward = findUpward(mesh, type, entities.data());
      PCU_ALWAYS_ASSERT(upward);
      return upward;
    }
  }

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

    edges.clear();
    oldTets.clear();
    oldVerts.clear();
    newVerts.clear();
    newLayer.clear();
    newEdges.clear();
    newTriangles.clear();
    elementsToChange.clear();
    old2New.clear();

    collectEdges(edge);
    computeElementSets();
    double oldMinQuality = collectOldTets();
    rebuildElements();

    bool valid = checkValidity(oldMinQuality);
    if (valid) destroyOldElements();
    else destroyNewElements();
    //TODO: check quality before accepting swap
    
    return valid;
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

  void LayerEdgeSwap::rebuildElements()
  {
    rebuildEdges();
    rebuildTriangles();
    rebuildTets();
    rebuildQuads();
    rebuildPrisms();
  }

  void LayerEdgeSwap::rebuildEdges()
  {
    for (size_t i=0; i<edges.size(); ++i)
      newEdges.push_back(mesh->createEntity(Mesh::EDGE, mesh->toModel(edges[i]), &newVerts[i*2]));
  }

  void LayerEdgeSwap::rebuildTriangles()
  {
    for (size_t e=0; e<edges.size(); ++e) {
      EntityArray triangles;
      mesh->getAdjacent(edges[e], Mesh::TRIANGLE, triangles);
      int v = 0;
      for (size_t t = 0; t < triangles.getSize(); ++t)
        if (getFlag(a, triangles[t], LAYER) && mesh->getType(triangles[t]) == Mesh::TRIANGLE) {
          Entity* oppEdge1 = getUpward(mesh, Mesh::EDGE, {oldVerts[e*2+v], newVerts[e*2]});
          Entity* oppEdge2 = getUpward(mesh, Mesh::EDGE, {oldVerts[e*2+v], newVerts[e*2+1]});
          Entity* triEdges[3]{newEdges[e], oppEdge1, oppEdge2};
          newTriangles.push_back(mesh->createEntity(Mesh::TRIANGLE, mesh->toModel(triangles[t]), triEdges));
          auto [it, inserted] = old2New.insert({triangles[t], newTriangles.back()});
          PCU_ALWAYS_ASSERT(inserted);
          v++;
        }
    }
  }

  void LayerEdgeSwap::rebuildTets() 
  {
    EntityVector oldFaces[2];
    for (size_t i = 0; i < oldTets.size(); ++i) {
      oldFaces[0].push_back(getTetFaceOppositeVert(mesh, oldTets[i], oldVerts[1]));
      oldFaces[1].push_back(getTetFaceOppositeVert(mesh, oldTets[i], oldVerts[0]));
    }
    PCU_ALWAYS_ASSERT(oldFaces[0].size() == 2);
    PCU_ALWAYS_ASSERT(oldFaces[1].size() == 2);

    Entity* oldFace = getTetFaceOppositeVert(mesh, oldTets[0], newVerts[0]);
    Entity* oppVert = getTriVertOppositeEdge(mesh, oldFace, edges[0]);
    Entity* edge1 = getUpward(mesh, Mesh::EDGE, {oppVert, newVerts[0]});
    Entity* edge2 = getUpward(mesh, Mesh::EDGE, {oppVert, newVerts[1]});
    Entity* triEdges[3] = {newEdges[0], edge1, edge2};
    Entity* newTetTriangle = mesh->createEntity(Mesh::TRIANGLE, mesh->toModel(oldFace), triEdges);
    for (size_t j=0; j<2; ++j) {
      Entity* tetFaces[4]{newTetTriangle, newTriangles[j], oldFaces[j][0], oldFaces[j][1]};
      newLayer.push_back(mesh->createEntity(Mesh::TET, mesh->toModel(oldTets[j]), tetFaces));
    }
  }

  void LayerEdgeSwap::rebuildQuads()
  {
    for (size_t i=0; i<edges.size()-1; ++i) {
      EntityArray quads;
      mesh->getAdjacent(edges[i], 2, quads);
      for (size_t q=0; q<quads.getSize(); q++) {
        if (mesh->getType(quads[q]) == Mesh::QUAD) {
          Entity* oppEdge1 = getUpward(mesh, Mesh::EDGE, {newVerts[i*2], newVerts[i*2+2]});
          Entity* oppEdge2 = getUpward(mesh, Mesh::EDGE, {newVerts[i*2+1], newVerts[i*2+3]});
          Entity* quadEdges[4]{newEdges[i], oppEdge1, newEdges[i+1], oppEdge2};
          Entity* newQuad = mesh->createEntity(Mesh::QUAD, mesh->toModel(quads[q]), quadEdges);
          auto [it, inserted] = old2New.insert({quads[q], newQuad});
          PCU_ALWAYS_ASSERT(inserted);
          break;
        }
      }
    }
  }

  void LayerEdgeSwap::rebuildPrisms()
  {
    for (size_t v=0; v<newVerts.size()-1; v+=2) {
      Entity* oppEdge1 = getUpward(mesh, Mesh::EDGE, {oldVerts[v], newVerts[v+1]});
      Entity* oppEdge2 = getUpward(mesh, Mesh::EDGE, {oldVerts[v+1], newVerts[v]});
      EntityArray quads1;
      mesh->getAdjacent(oppEdge1, 2, quads1);
      EntityArray quads2;
      mesh->getAdjacent(oppEdge2, 2, quads2);
      Entity* oppQuad1;
      Entity* oppQuad2;
      for (size_t q=0; q<quads1.size(); ++q)
      if (mesh->getType(quads1[q]) == Mesh::QUAD && old2New.find(quads1[q]) == old2New.end())
        oppQuad1 = quads1[q];
      for (size_t q=0; q<quads2.size(); ++q)
      if (mesh->getType(quads2[q]) == Mesh::QUAD && old2New.find(quads2[q]) == old2New.end())
        oppQuad2 = quads2[q];
      old2New[oppQuad1] = oppQuad2;
      old2New[oppQuad2] = oppQuad1;
    }

    APF_ITERATE(EntitySet, elementsToChange, it) {
      if (mesh->getType(*it) == Mesh::TET) continue;
      Entity* prismFaces[5];
      mesh->getDownward(*it, 2, prismFaces);
      for (int i=0; i<5; i++) {
        auto newFaceIt = old2New.find(prismFaces[i]);
        if (newFaceIt != old2New.end()) prismFaces[i] = newFaceIt->second;
      }
      newLayer.push_back(mesh->createEntity(Mesh::PRISM, mesh->toModel(*it), prismFaces));
      setFlagOnClosure(a, newLayer.back(), LAYER);
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