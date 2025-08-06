#include "maAdapt.h"

namespace ma {

  enum FaceSwapType {
    Two2Two,
    Two2Three,
  };

  class FaceSwap
  {
    public:
    FaceSwap(Adapt* a, Entity* f):
      adapt(a),
      mesh(a->mesh),
      face(f)
    {}

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
      apf::Up adjTets;
      mesh->getUp(face, adjTets);
      Entity* newEdgeVerts[2];
      newEdgeVerts[0] = getTetVertOppositeTri(mesh, adjTets.e[0], face);
      newEdgeVerts[1] = getTetVertOppositeTri(mesh, adjTets.e[1], face);

      //TODO: determine two2two or two2three by checking if two adjacent sides are coplanar
      //        check that they are coplanar by comparing normals
      //TODO: check quality from swap depending each case

      if (type == Two2Two) {

      }
      else if (type == Two2Three) {
        
      }
      
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
  };

  bool runEdgeSwap(Adapt* a, Entity* face)
  {
    FaceSwap faceSwap(a, face);
    if (!faceSwap.topoCheck())
      return false;
    if (!faceSwap.geomCheck())
      return false;
    if (!faceSwap.sizeCheck())
      return false;
    return true;
  }

}