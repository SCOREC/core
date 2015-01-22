#ifndef PARMA_DCPART_H_
#define PARMA_DCPART_H_

#include "apf.h"
#include "apfMesh.h"
#include <vector>
#include <map>

class dcPart {
   public:
      dcPart(apf::Mesh*& mesh, unsigned verbose=0);
      ~dcPart();
      unsigned numDisconnectedComps();
      bool isIsolated(apf::MeshEntity* e);
      void fix();
   private:
      dcPart() {}
      unsigned walkPart(unsigned visited);
      void markIsolated(const unsigned dcComp);
      unsigned maxContactNeighbor(const unsigned dcComp);
      void setupPlan(std::map<unsigned,unsigned> & dcCompTgts, 
          apf::Migration* plan);
      int totNumDc();

      std::vector<unsigned> dcCompSz;
      std::vector<unsigned> dcCompNbor;
      apf::MeshTag* vtag;
      apf::MeshTag* isotag;
      apf::Mesh* m;
      unsigned verbose;
};

#endif
