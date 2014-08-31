#ifndef PARMA_PARTINFO_H_
#define PARMA_PARTINFO_H_

#include "apf.h"
#include "apfMesh.h"
#include "parma_priority.h"
#include <vector>
#include <set>

class partInfo {
   public: 
      int id;
      double weight[4];
      std::vector<int> isCandidateAp;       
      std::vector<int> adjPartIds;
      std::vector< apf::Array<double, 4> > adjPartWeights;
      bool isCandidate(const int adjPid);
      void print(int key);
      partInfo(apf::Mesh* m, apf::MeshTag* wtag, const priorityList& pl, const int plIdx, double* maxImbW);
      void getGreedyMigrationSchedule2(const priorityList &pl, const int plIdx, double* maxImbW);
   private:
      partInfo();
      void get(apf::Mesh* m, apf::MeshTag* wtag);
      void getGreedyMigrationSchedule(const priorityList &pl, const int plIdx, double* maxImbW);
      void initNeighbors(std::set<int>& ap);
      int getAdjPartWeight(apf::Mesh* m);
      int sendWeightToNeighbors();
      int recvWeightFromNeighbors();
};

#endif 
