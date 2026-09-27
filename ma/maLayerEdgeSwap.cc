#include "maAdapt.h"
#include "maLayerEdgeSwap.h"
#include "maLayer.h"
#include "maCrawler.h"
#include "pcu_util.h"
#include "maDBG.h"
#include "maFirstProblemPlane.h"

namespace ma {

  LayerEdgeSwap::LayerEdgeSwap(Adapt* adapt) {
    a = adapt;
    mesh = a->mesh;
  }

  bool LayerEdgeSwap::run(Entity* edge)
  {
    if (!a->hasLayer) return false;
    if (!a->input->shouldSwapLayer) return false;
    if (isOnModelEdge(mesh, edge)) return false;
    if (!getFlag(a, edge, LAYER)) return false;
    //TODO: only accept edge if its the layer base

    edges.clear();
    oldVerts.clear();
    newVerts.clear();
    newLayer.clear();
    newEdges.clear();
    newTriangles.clear();
    elementsToChange.clear();

    collectEdges(edge);
    computeElementSets();
    rebuildElements();
    printInfo();
    destroyOldElements();

    //TODO: check quality before accepting swap
    
    for (int i=0; i<edges.size(); i++) clearFlag(a, edges[i], CHECKED);
    return true;
  }

  void LayerEdgeSwap::collectEdges(Entity* edge) {
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
          oldTriangles.push_back(faces[i]);
          newVerts.push_back(getTriVertOppositeEdge(mesh,faces[i],edge));
        }
      
      edge = getOtherEdge(mesh, edge, p);
      if (!edge) break;
    }
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

  void LayerEdgeSwap::rebuildTets() 
  {
    //TODO: pass in old tets, new edge and triangles instead of global
    EntityVector oldFaces[2];
    EntityArray tets;
    mesh->getAdjacent(edges[0], 3, tets);
    for (size_t i=0; i<tets.getSize(); ++i)
    if (mesh->getType(tets[i]) == Mesh::TET) {
      oldFaces[0].push_back(getTetFaceOppositeVert(mesh, tets[i], oldVerts[1]));
      oldFaces[1].push_back(getTetFaceOppositeVert(mesh, tets[i], oldVerts[0]));
    }
    PCU_ALWAYS_ASSERT(oldFaces[0].size() == 2);
    PCU_ALWAYS_ASSERT(oldFaces[1].size() == 2);

    Entity* oldFace = getTetFaceOppositeVert(mesh, tets[0], newVerts[0]);
    Entity* oppVert = getTriVertOppositeEdge(mesh, oldFace, edges[0]);
    Entity* edgeVerts1[2] = {oppVert, newVerts[0]};
    Entity* edgeVerts2[2] = {oppVert, newVerts[1]};
    Entity* edge1 = findUpward(mesh, Mesh::EDGE, edgeVerts1);
    Entity* edge2 = findUpward(mesh, Mesh::EDGE, edgeVerts2);
    if (!edge1 || !edge2) printf("TRI EDGE NOT FOUND\n");
    Entity* triEdges[3] = {newEdges[0], edge1, edge2};
    Entity* newVerticalTriangle = mesh->createEntity(mesh->getType(oldFace), mesh->toModel(oldFace), triEdges);
    for (size_t j=0; j<2; ++j) {
      Entity* tetFaces[4]{newVerticalTriangle, newTriangles[j], oldFaces[j][0], oldFaces[j][1]};
      newLayer.push_back(mesh->createEntity(Mesh::TET, mesh->toModel(tets[j]), tetFaces));
    }
  }

  void LayerEdgeSwap::rebuildElements()
  {
    EntityVector newQuads;
    EntityVector oldQuads;
    EntityVector keptEdges;
    EntityVector keptQuads;
    for (size_t i=0; i<edges.size(); ++i){
      newEdges.push_back(mesh->createEntity(mesh->getType(edges[i]), mesh->toModel(edges[i]), &newVerts[i*2]));

      for (size_t j = 0; j < 2; ++j) {
        Entity* oppEdges[2];
        Entity* av1[2] = {oldVerts[(i*2)+j], newVerts[i*2]};
        oppEdges[0] = findUpward(mesh, Mesh::EDGE, av1);
        Entity* av2[2] = {oldVerts[(i*2)+j], newVerts[(i*2)+1]};
        oppEdges[1] = findUpward(mesh, Mesh::EDGE, av2);
        if (!oppEdges[0] || !oppEdges[1]) printf("EDGE NOT FOUND\n");
        Entity* triEdges[3]{newEdges[i], oppEdges[0], oppEdges[1]};
        newTriangles.push_back(mesh->createEntity(mesh->getType(oldTriangles[i]), mesh->toModel(oldTriangles[i]), triEdges));

        if (i%2 != 0) continue;
        for (size_t k = 0; k < 2; ++k) {
          EntityArray faces;
          mesh->getAdjacent(oppEdges[k], 2, faces);
          for (size_t l = 0; l < faces.getSize(); ++l)
            if (mesh->getType(faces[l]) == Mesh::QUAD){
              keptQuads.push_back(faces[l]);
            }
        }
      }

      if (i%2 != 0) continue;
      EntityArray faces;
      mesh->getAdjacent(edges[i], 2, faces);
      for (size_t j = 0; j < faces.getSize(); ++j)
        if (mesh->getType(faces[j]) == Mesh::QUAD)
          oldQuads.push_back(faces[j]);
    }

    for (size_t i=0; i<newVerts.size()-2; ++i){
      Entity* av[2] = {newVerts[i], newVerts[i+2]};
      Entity* keptEdge = findUpward(mesh, Mesh::EDGE, av);
      if (!keptEdge) printf("EDGE NOT FOUND\n");
      else keptEdges.push_back(keptEdge);
    }

    for (size_t i=0; i<oldQuads.size(); ++i) {
      Entity* faceEdges[4]{newEdges[i], keptEdges[i*2], newEdges[i+1], keptEdges[i*2+1]};
      newQuads.push_back(mesh->createEntity(mesh->getType(oldQuads[i]), mesh->toModel(oldQuads[i]), faceEdges));
    }

    rebuildTets();

    // int kq = 0;
    // int nq = 0;
    // for (size_t i=0; i<newQuads.size(); ++i) {
    //   EntityArray prisms;
    //   mesh->getAdjacent(oldQuads[i], 3, prisms);
    //   for (size_t j=0; j<prisms.getSize(); ++j) {
    //     Entity* prismFaces[5]{newTriangles[nq], newTriangles[nq+2], newQuads[i], keptQuads[kq+1], keptQuads[kq]};
    //     newLayer.push_back(mesh->createEntity(mesh->getType(prisms[j]), mesh->toModel(prisms[j]), prismFaces));
    //     nq++;
    //     kq+=2;
    //   }
    // }

    ma_dbg::useFieldInfo(a, [this, keptEdges, keptQuads, newQuads] {
      EntityVector keptE = keptEdges;
      EntityVector copyQuads = keptQuads;
      EntityVector nQuads = newQuads;
      ma_dbg::flagEntity(a, 1, "edges_to_swap", &edges[0], edges.size());
      ma_dbg::flagEntity(a, 0, "target_verts", &newVerts[0], newVerts.size());
      ma_dbg::flagEntity(a, 1, "kept_edges", &keptE[0], 2);
      ma_dbg::flagEntity(a, 1, "new_edges", &newEdges[0], 2);
      ma_dbg::flagEntity(a, 2, "kept_quads", &copyQuads[0], 2);
      ma_dbg::flagEntity(a, 2, "new_triangles", &newTriangles[0]+2, 1);
      ma_dbg::flagEntity(a, 2, "new_quads", &nQuads[0], 1);
      ma_dbg::flagEntity(a, 3, "new_layer", &newLayer[0], newLayer.size());

      apf::writeVtkFiles("mesh_tets", mesh, 3);
      apf::writeVtkFiles("mesh_faces", mesh, 2);
      apf::writeVtkFiles("mesh_edges", mesh, 1);
    });
    printf("COLLECTED\n");
    exit(1);
  }

  void LayerEdgeSwap::destroyOldElements()
  {
    APF_ITERATE(EntitySet, elementsToChange, it)
      destroyElement(a, *it);
  }
  
  void LayerEdgeSwap::printInfo() {
    ma_dbg::useFieldInfo(a, [this] {
      ma_dbg::flagEntity(a, 1, "edges_to_swap", &edges[0], edges.size());
      ma_dbg::flagEntity(a, 0, "target_verts", &newVerts[0], newVerts.size());

      apf::writeVtkFiles("mesh_tets", mesh, 3);
      apf::writeVtkFiles("mesh_faces", mesh, 2);
      apf::writeVtkFiles("mesh_edges", mesh, 1);
    });
    printf("COLLECTED\n");
    exit(1);
  }
}