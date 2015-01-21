#ifndef PARMA_DCPART_H_
#define PARMA_DCPART_H_

#include "apf.h"
#include "apfMesh.h"
#include <vector>
#include <set>
#include <map>

// < componentIdx, mergeTgtPartId >
typedef std::map<size_t,int> migrTgt;

class dcPart {
   public:
      dcPart(apf::Mesh*& mesh, unsigned verbose=0);
      ~dcPart();
      int numDisconnectedComps();
      void makeDisconnectedComps(const int numDcComps);
      bool isIsolated(apf::MeshEntity* e);
      void fix();
   private:
      dcPart() {}
      size_t walkPart(size_t visited);
      void markIsolated(const size_t dcComp);
      unsigned maxContactNeighbor(const size_t dcComp);
      void setupPlan(migrTgt& dcCompTgts, apf::Migration* plan);
      int totNumDc();

      std::vector<size_t> dcCompSz;
      std::vector<unsigned> dcCompNbor;
      apf::MeshTag* vtag;
      apf::MeshTag* isotag;
      apf::Mesh* m;
      unsigned verbose;
};

#endif
