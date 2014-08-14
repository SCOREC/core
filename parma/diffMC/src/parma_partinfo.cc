#include "parma_partinfo.h"
#include "parma_commons.h"
#include "PCU.h"
#include <stdio.h>
#include <sstream>

using std::set;

using apf::Mesh;
using apf::MeshEntity;
using apf::MeshIterator;
using apf::MeshTag;

using parmaCommons::isLess;
using parmaCommons::isMore;
using parmaCommons::status;

int partInfo::sendWeightToNeighbors() {
  APF_ITERATE(std::vector<int>,adjPartIds,adjPartIdItr) {
    const int destRank = *adjPartIdItr;  
    //pack weight for current part
    PCU_COMM_PACK(destRank, weight);
  }
  return 0;
}

int partInfo::recvWeightFromNeighbors() {
  while (PCU_Comm_Listen()) {
    int srcRank;
    PCU_Comm_From(&srcRank);
    int found = 0;
    for (size_t apIdx = 0; apIdx < adjPartIds.size(); apIdx++) {
      if (adjPartIds[apIdx] == srcRank) {
        found = 1;
        for(int entDim=0; entDim<4; entDim++) {
          double weight;
          PCU_COMM_UNPACK(weight);
          adjPartWeights[apIdx][entDim] = weight;
        }
      }
    } 
    assert(1 == found);
  }
  return 0;
}

/**
 * @brief for each local part get the id of face-adjacent parts and the weight they have 
 * @param mesh (In) partitioned mesh
 * @param part (InOut) partInfo object
 * @return zero on success, non-zero otherwise
 */
int partInfo::getAdjPartWeight(Mesh* m) {
    if ( ! PCU_Comm_Initialized() ) 
       PCU_Comm_Init();
    PCU_Comm_Start(PCU_GLOBAL_METHOD);

    const int dim = m->getDimension() - 1;
    assert(dim == 2);
    set<int> adjParts;
    getFacePeers(m, adjParts);
    initNeighbors(adjParts);

    int ierr = sendWeightToNeighbors();
    if (ierr != 0) return ierr;
    PCU_Comm_Send();
    ierr = recvWeightFromNeighbors();
    if (ierr != 0) return ierr;

    return 0;
}

inline void getWeight(apf::Mesh* m, apf::MeshTag* wtag, double* w) {
   apf::MeshEntity* e;
   for(int d=0; d<=m->getDimension(); d++) { 
      w[d] = 0;
      MeshIterator* itr = m->begin(d);
      while( (e = m->iterate(itr)) ) {
         double weight = 1;
         if ( m->hasTag(e, wtag) ) {
           m->getDoubleTag(e, wtag, &weight);
         }
	 w[d] += weight;
      }
      m->end(itr);
   }
}

/**
 * @brief determine if a part is a candidate for migration
 * @remark if the adjacent part has weight greater than 97% of the source 
 *         parts weight, for any entity type of equal or higher priority 
 *         then it is not a candidate for migration
 *
 * @param pl (In) priority list 
 * @param plIdx (In) priority list idx
 */
void partInfo::getGreedyMigrationSchedule(const priorityList &pl, const int plIdx, double* maxImbW ) {
   for ( size_t apIdx=0; apIdx<adjPartIds.size(); apIdx++) {
      int ic = 1;
      for ( int pIdx=plIdx; pIdx>=0; pIdx-- ) {
	 const int eDim = pl.entDim[pIdx];
	 if ( isLess(weight[eDim], maxImbW[eDim])
	       || isMore(adjPartWeights[apIdx][eDim], 0.97*weight[eDim]) ) {
	    ic = 0;
	 }
      }
      isCandidateAp[apIdx] = ic;
   }
}

void partInfo::getGreedyMigrationSchedule2(const priorityList &pl, const int plIdx, double* maxImbW ) {
   for ( size_t apIdx=0; apIdx<adjPartIds.size(); apIdx++) {
      int ic = 1;
      const int eDim = pl.entDim[plIdx];
      if ( isLess(weight[eDim], maxImbW[eDim])
          || isMore(adjPartWeights[apIdx][eDim], 0.97*weight[eDim]) ) {
        ic = 0;
      }
      isCandidateAp[apIdx] = ic;
   }
}

void partInfo::get(Mesh* m, MeshTag* wtag) {
   id = m->getId();
   double w[4] = {0,0,0,0};
   getWeight(m, wtag, w);
   for(int i=0; i<4; i++)
      weight[i] = w[i];

   int ierr = getAdjPartWeight(m);
   if (ierr != 0) exit(1);
}

void partInfo::initNeighbors(set<int>& ap) {
  const int numNp = ap.size();
  assert(numNp > 0);
/* the clang static analyzer emits a warning in GNU STL
   code when resizing vectors of non-trivial things.
   this swap hack does the same without a warning */
//adjPartWeights.resize(numNp);
  std::vector< apf::Array<double, 4> > clang_v(numNp);
  std::swap(adjPartWeights, clang_v);
  APF_ITERATE(set<int>, ap, it) {
    adjPartIds.push_back(*it);
  }
  isCandidateAp.resize(numNp);
  for(int i=0; i<numNp; i++) {
    isCandidateAp[i] = 0;
    for(int j=0; j<4; j++) 
      adjPartWeights[i][j] = 0;
  }
}

bool partInfo::isCandidate(const int adjPid) {
   for(size_t i=0; i<adjPartIds.size(); i++) {
      if( adjPid == adjPartIds[i] && 1 == isCandidateAp[i] ) 
         return true; 
   }
   return false;
}

void partInfo::print(int key) {
   std::stringstream msg;
   msg << key << " id weight " << id << ' '
       << weight[0] << ' ' 
       << weight[1] << ' ' 
       << weight[2] << ' ' 
       << weight[3] << '\n';
   for(size_t i=0; i<adjPartIds.size(); i++) {
      msg << key << " id apId isCandidate weight " << id << ' '
            << adjPartIds[i] << ' '
            << isCandidateAp[i] << ' '
            << adjPartWeights[i][0] << ' '
            << adjPartWeights[i][1] << ' '
            << adjPartWeights[i][2] << ' '
            << adjPartWeights[i][3] << '\n';
   }
   std::string s = msg.str();
   status(s.c_str());
}

partInfo::partInfo(Mesh* m, MeshTag* wtag, const priorityList& pl, const int plIdx, double* maxImb) {
   id = -1;
   for(int i=0; i<4; i++) 
      weight[i] = -1;
   get(m, wtag);
   getGreedyMigrationSchedule(pl, plIdx, maxImb);
}

