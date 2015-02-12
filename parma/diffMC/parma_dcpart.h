#ifndef PARMA_DCPART_H_
#define PARMA_DCPART_H_

#include "apf.h"
#include "apfMesh.h"
#include <vector>

class dcPart {
   public:
      dcPart(apf::Mesh*& mesh, unsigned verbose=0);
      ~dcPart();
      unsigned getNumComps();
      unsigned getNumIso();
      unsigned getCompSize(unsigned i);
      unsigned getCompPeer(unsigned i);
      unsigned numDisconnectedComps();
      bool isIsolated(apf::MeshEntity* e);
      unsigned compId(apf::MeshEntity* e);
      apf::MeshEntity* getSeedEnt(unsigned i);
      void fix();
   protected:
      void reset();
   private:
      dcPart() {}
      unsigned walkPart(unsigned visited);
      void markIsolated(const unsigned dcComp);
      unsigned maxContactNeighbor(const unsigned dcComp);

      unsigned numIso;
      std::vector<unsigned> dcCompSz;
      std::vector<unsigned> dcCompNbor;
      apf::MeshTag* vtag;
      apf::MeshTag* isotag;
      apf::Mesh* m;
      unsigned verbose;
};

class dcPartFixer { 
  public:
    dcPartFixer(apf::Mesh* mesh, unsigned verbose=0);
    ~dcPartFixer();
  private:
    dcPartFixer();
    class PartFixer;
    PartFixer* pf;
};

namespace parma {
  class dcComponents {
    public:
      dcComponents(apf::Mesh* mesh, unsigned verbose=0);
      ~dcComponents();

      unsigned size();
      unsigned numIso();

      bool has(apf::MeshEntity* e);
      unsigned getId(apf::MeshEntity* e);

      apf::MeshEntity* getCore(unsigned i);

      bool bdryHas(unsigned i, apf::MeshEntity* e);
      void beginBdry(unsigned i);
      apf::MeshEntity* iterateBdry();
      void endBdry();
    private:
      class Components;
      Components* c;
      class BdryItr;
      BdryItr* bItr;
  };
}

#endif
