#ifndef PARMA_COMPONENTS_H_
#define PARMA_COMPONENTS_H_

#include "parma_dcpart.h"

namespace parma {
  typedef std::set<apf::MeshEntity*> Level;
  class dcComponents::Components : public dcPart {
    public:
      Components(apf::Mesh* mesh, unsigned verbose=0);
      ~Components();

      unsigned size();
      unsigned iso();

      unsigned getDepth(unsigned i);
      void setDepth(unsigned i, unsigned d);

      Level* getBdry(unsigned i);
      Level* getCore(unsigned i);
      apf::MeshEntity* getCoreVtx(unsigned i);

      bool has(apf::MeshEntity* e);
      unsigned getId(apf::MeshEntity* e);
    private:
      Components();
      void markVertices();
      void walkComp(apf::MeshEntity* src, unsigned comp);
      void setElmVtxIds(apf::Downward& verts, const int nv, unsigned comp);
      void addElmVtxToBdry(apf::Downward& verts, const int nv, unsigned comp);
      void setId(apf::MeshEntity* e, unsigned id);
      void walkInward(unsigned compId);
      void getCoreVerts();
      void getCoreVtx();
      void sortByDepth();
      void reorder(unsigned*);
      apf::Mesh* m;
      unsigned vb;
      apf::MeshTag* idT;
      unsigned n;
      unsigned* id;
      unsigned* depth;
      Level* bdry;
      Level* core;
  };
}

#endif
