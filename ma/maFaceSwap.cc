#include "maAdapt.h"
#include "maSnapper.h"
#include "apfGeometry.h"
#include "maShapeHandler.h"
#include "maDBG.h"

//TODO: add comments

namespace ma {

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
    FaceSwap(Adapt* a, Entity* f, bool improve) : 
      adapt(a), mesh(a->mesh), face(f), improveQuality(improve)
    {
      cavity.init(a);
      numNewTets=0;
    }

    void printCavityBefore()
    {
      EntityArray cavity;
      for (int i=0; i<2; i++){
        cavity.append(oldTets[i]);
      }
      ma_dbg::createCavityMesh(adapt, cavity, "BeforeFaceSwap");
    }

    void printCavityAfter()
    {
      EntityArray cavity;
      for (int i=0; i<numNewTets; i++){
        cavity.append(newTets[i]);
      }
      ma_dbg::createCavityMesh(adapt, cavity, "AfterFaceSwap");
    }

    bool invalid(Entity* verts[4])
    {
      Vector v[4] {getPosition(mesh, verts[0]), getPosition(mesh, verts[1]), 
                  getPosition(mesh, verts[2]), getPosition(mesh, verts[3])};
      if ((apf::cross((v[1] - v[0]), (v[2] - v[0])) * (v[3] - v[0])) < 0)
        return true;
      return false;
    }

    bool invalidSwap() {
      Entity* opp1 = getTetVertOppositeTri(mesh, oldTets[0], face);
      Entity* opp2 = getTetVertOppositeTri(mesh, oldTets[1], face);
      Entity* faceVerts[3];
      mesh->getDownward(face, 0, faceVerts);
      Entity* tetVertsA[4] = {faceVerts[0], faceVerts[1], opp1, opp2};
      if (invalid(tetVertsA)) return true;
      Entity* tetVertsB[4] = {faceVerts[1], faceVerts[2], opp1, opp2};
      if (invalid(tetVertsB)) return true;
      Entity* tetVertsC[4] = {faceVerts[2], faceVerts[0], opp1, opp2};
      if (invalid(tetVertsC)) return true;
      return false;
    }

    bool topoCheck()
    {   
      int modelDim = mesh->getModelType(mesh->toModel(face));
      if (modelDim == 2)
        return false;
      if (getFlag(adapt,face,DONT_SWAP))
        return false;
      mesh->getAdjacent(face, 3, oldTets);
      if (invalidSwap())
        return false;
      return true;
    }

    void findNumNewTets()
    {
      numNewTets = 3;
      Entity* faceVerts[3];
      mesh->getDownward(face, 0, faceVerts);
      for (int i=0; i<3; i++) {
        Entity* face0 = getTetFaceOppositeVert(mesh, oldTets[0], faceVerts[i]);
        Entity* face1 = getTetFaceOppositeVert(mesh, oldTets[1], faceVerts[i]);
        Vector normal0 = getTriNormal(mesh, face0);
        Vector normal1 = getTriNormal(mesh, face1);
        if (apf::areClose(normal0, normal1, 1e-10)) {
          numNewTets = 0; //TODO: TEST Two2Two CASE
          commonEdge = findCommonEdge(mesh, face0, face1);
        }
      }
    }

    void buildTwo2Two()
    {
      printf("\t[Warning]: found two2two face swap, this is unteasted please test\n");
      Model* model = mesh->toModel(oldTets[0]);
      Entity* commEdgeVerts[2];
      mesh->getDownward(commonEdge, 0, commEdgeVerts);

      Entity* newTetVerts0[4];
      Entity* newTetVerts1[4];
      newTetVerts0[0] = commEdgeVerts[0];
      newTetVerts1[0] = commEdgeVerts[1];

      fillAdjVerts(mesh, newTetVerts0, oldTets);
      fillAdjVerts(mesh, newTetVerts1, oldTets);

      newTets[0] = buildElement(adapt, model, apf::Mesh::TET, newTetVerts0);
      newTets[1] = buildElement(adapt, model, apf::Mesh::TET, newTetVerts1);
    }

    void buildTwo2Three()
    {
      Model* model = mesh->toModel(oldTets[0]);
      Entity* newTetVerts[4];
      newTetVerts[0] = getTetVertOppositeTri(mesh, oldTets[0], face);
      newTetVerts[1] = getTetVertOppositeTri(mesh, oldTets[1], face);

      Entity* faceEdges[3];
      mesh->getDownward(face, 1, faceEdges);
      for (int i=0; i<3; i++) {
        Entity* faceEdgeVerts[2];
        mesh->getDownward(faceEdges[i], 0, faceEdgeVerts);
        newTetVerts[2]=faceEdgeVerts[0];
        newTetVerts[3]=faceEdgeVerts[1];
        newTets[i] = buildElement(adapt, model, apf::Mesh::TET, newTetVerts);
      }
    }

    void destroyOldElements()
    {
      cavity.transfer(oldTets);
      destroyElement(adapt, oldTets[0]);
      destroyElement(adapt, oldTets[1]);
    }

    void cancel()
    {
      for (int i=0; i<numNewTets; i++)
        destroyElement(adapt, newTets[i]);
    }

    bool geomCheck()
    {
      findNumNewTets();
      // printCavityBefore();
      cavity.beforeBuilding();
      if (numNewTets == 0)
        return false;
      else if (numNewTets == 2)
        buildTwo2Two();
      else if (numNewTets == 3)
        buildTwo2Three();
      cavity.afterBuilding();
      cavity.fit(oldTets);
      // printCavityAfter();
      return true;
    }

    bool sizeCheck()
    {
      if (!improveQuality) return true;
      for (int i=0; i<numNewTets; i++) {
        Entity* edges[6];
        mesh->getDownward(newTets[i], 1, edges);
        for (int e=0; e<6; e++)
          if (adapt->sizeField->measure(edges[e]) > MAXLENGTH)
            return false;
      }
      return true;
    }

    bool qualityCheck()
    {
      double qualityToBeat = adapt->input->validQuality;
      if (improveQuality)
        qualityToBeat = std::max(getWorstQuality(adapt, oldTets), adapt->input->validQuality);

      for (int i=0; i<numNewTets; i++) {
        double quality = adapt->shape->getQuality(newTets[i]);
        if (quality < qualityToBeat)
          return false;
      }
      return true;
    }

    private:
    Adapt* adapt;
    Mesh* mesh;
    Entity* face;
    Cavity cavity;
    int numNewTets;
    Upward oldTets;
    Entity* newTets[3];
    bool improveQuality;
    Entity* commonEdge;
  };

  bool runFaceSwap(Adapt* a, Entity* face, bool improveQuality)
  {
    FaceSwap faceSwap(a, face, improveQuality);

    if (faceSwap.topoCheck() 
      && faceSwap.geomCheck() 
      && faceSwap.sizeCheck() 
      && faceSwap.qualityCheck()) {
      faceSwap.destroyOldElements();
      return true;
    }
    else {
      faceSwap.cancel();
      return false;
    }
  }

}