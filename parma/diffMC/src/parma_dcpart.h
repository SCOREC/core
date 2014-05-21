#ifndef PARMA_DCPART_H_
#define PARMA_DCPART_H_

#include "apf.h"
#include "apfMesh.h"
#include <vector>
#include <set>
#include <map>

typedef std::map<int,int> migrTgt;

class dcPart {
   public:
      dcPart(apf::Mesh*& mesh);
      ~dcPart();
      int numDisconnectedComps();
      void makeDisconnectedComps(const int numDcComps);
      void fix();
   private:
      dcPart();
      void init(apf::Mesh*& mesh);
      int walkPart(int visited);
      int checkResidence(const int dcComp);
      void setupPlan(migrTgt& dcCompTgts, apf::Migration* plan);

      std::vector<int> dcCompSz;
      apf::MeshTag* vtag;
      apf::Mesh* m;
};

#endif
