#include "maAdapt.h"
#include "maSnapper.h"
#include "apfGeometry.h"
#include "maDBG.h"

//TODO: move useful functions to header
//TODO: break up functions into smaller chunks
//TODO: add comments

namespace ma {

  enum FaceSwapType {
    Two2Two,
    Two2Three,
  };

  Entity* findCommonEdge(Mesh* mesh, Entity* face1, Entity* face2)
  {
    Entity* face1Edges[3];
    mesh->getDownward(face1, 1, face1Edges);
    for (int i=0; i<3; i++)
      if (isLowInHigh(mesh, face2, face1Edges[i]))
        return face1Edges[i];
    return face1;
  }

  void fillAdjVerts(Mesh* mesh, Entity* newTetVerts[4], const Upward& adjTets)
  {
    apf::Up adjEdges;
    mesh->getUp(newTetVerts[0], adjEdges);
    int j=1;
    for (int i=0; i<adjEdges.n; i++) {
      if (isLowInHigh(mesh, adjTets[0] , adjEdges.e[i]))
        newTetVerts[j++] = apf::getEdgeVertOppositeVert(mesh, adjEdges.e[i], newTetVerts[0]);
    }
  }

  class FaceSwap
  {
    public:
    FaceSwap(Adapt* a, Entity* f): adapt(a), mesh(a->mesh), face(f)
    {
      cavity.init(a);
    }

    bool topoCheck()
    {   
      int modelDim = mesh->getModelType(mesh->toModel(face));
      if (modelDim == 2)
        return false;
      if (getFlag(adapt,face,DONT_SWAP))
        return false;
      return true;
    }

    bool geomCheck()
    {
      Upward adjTets;
      mesh->getAdjacent(face, 3, adjTets);

      type = Two2Three;
      Entity* faceVerts[3];
      mesh->getDownward(face, 0, faceVerts);
      Entity* commonEdge;
      for (int i=0; i<3; i++) {
        Entity* face0 = getTetFaceOppositeVert(mesh, adjTets[0], faceVerts[i]);
        Entity* face1 = getTetFaceOppositeVert(mesh, adjTets[1], faceVerts[i]);
        Vector normal0 = getTriNormal(mesh, face0);
        Vector normal1 = getTriNormal(mesh, face1);
        if (apf::areClose(normal0, normal1, 1e-10)) {
          type = Two2Two;
          commonEdge = findCommonEdge(mesh, face0, face1);
        }
      }

      //TODO: check quality from swap depending each case
      cavity.beforeBuilding();
      Model* model = mesh->toModel(adjTets[0]);
      if (type == Two2Two) {
        printf("\t[Warning]: found two2two face swap, this is unteasted please test\n");
        Entity* commEdgeVerts[2];
        mesh->getDownward(commonEdge, 0, commEdgeVerts);

        Entity* newTetVerts0[4];
        Entity* newTetVerts1[4];
        newTetVerts0[0] = commEdgeVerts[0];
        newTetVerts1[0] = commEdgeVerts[1];

        fillAdjVerts(mesh, newTetVerts0, adjTets);
        fillAdjVerts(mesh, newTetVerts1, adjTets);

        buildElement(adapt, model, apf::Mesh::TET, newTetVerts0);
        buildElement(adapt, model, apf::Mesh::TET, newTetVerts1);
      }
      else if (type == Two2Three) {
        Entity* newTetVerts[4];
        newTetVerts[0] = getTetVertOppositeTri(mesh, adjTets[0], face);
        newTetVerts[1] = getTetVertOppositeTri(mesh, adjTets[1], face);

        Entity* faceEdges[3];
        mesh->getDownward(face, 1, faceEdges);
        for (int i=0; i<3; i++) {
          Entity* faceEdgeVerts[2];
          mesh->getDownward(faceEdges[i], 0, faceEdgeVerts);
          newTetVerts[2]=faceEdgeVerts[0];
          newTetVerts[3]=faceEdgeVerts[1];
          buildElement(adapt, model, apf::Mesh::TET, newTetVerts);
        }
      }
      cavity.afterBuilding();
      cavity.fit(adjTets);
      cavity.transfer(adjTets);
      destroyElement(adapt, adjTets[0]);
      destroyElement(adapt, adjTets[1]);
      return true;
    }

    bool sizeCheck()
    {
      if (type == Two2Two) {

      }
      else if (type == Two2Three) {

      }
      return true;
    }

    private:
    Mesh* mesh;
    Adapt* adapt;
    Entity* face;
    FaceSwapType type;
    Cavity cavity;
  };

  bool runFaceSwap(Adapt* a, Entity* face)
  {
    printf("InFaceSwap\n");
    FaceSwap faceSwap(a, face);
    if (!faceSwap.topoCheck())
      return false;
    printf("TopoCheck\n");
    if (!faceSwap.geomCheck())
      return false;
    printf("GeomCheck\n");
    exit(1);
    if (!faceSwap.sizeCheck())
      return false;
    return true;
  }

}