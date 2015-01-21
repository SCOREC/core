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
      void fix();
   private:
      dcPart() {}
      void init(apf::Mesh*& mesh);
      size_t walkPart(size_t visited);
      unsigned maxContactNeighbor(const size_t dcComp);
      void setupPlan(migrTgt& dcCompTgts, apf::Migration* plan);
      int totNumDc();

      std::vector<size_t> dcCompSz;
      std::vector<unsigned> dcCompNbor;
      apf::MeshTag* vtag;
      apf::Mesh* m;
      unsigned verbose;
};

#endif
