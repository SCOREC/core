#include "PCU.h"
#include "parma_diffmc.h"
#include "parma_commons.h"
#include "parma_imbinfo.h"
#include "parma_dcpart.h"
#include "parma_hist.h"
#include <stdio.h>
#include <vector>
#include <set>
#include <list>
#include <sstream>

#include "mpi.h"

#include "parma_meshaux.h"

using namespace apf;

using parmaCommons::debug;
using parmaCommons::status;
using parmaCommons::error;
using parmaCommons::isLess;
using parmaCommons::isEql;
using parmaCommons::printElapsedTime;

using std::vector;
using std::set;
using std::list;

namespace {

inline void markAsVisited(Mesh* m, MeshTag* vtag, int visited, MeshEntity* rgn, Downward& adjFaces, const int numFaces) {
   for(int i=0; i<numFaces; i++)
      m->setIntTag(adjFaces[i], vtag, &visited); 
   m->setIntTag(rgn, vtag, &visited); 
}

inline double addToPlan(Mesh* m, MeshTag* wtag, 
      MeshEntity* rgn, const int destPid, Migration* plan) {
   plan->send(rgn, destPid);
   double elmW = 0;
   m->getDoubleTag(rgn, wtag, &elmW);
   return elmW;
}

/**
 * @brief determine if the elements in rgns are dim-1 (face in 3D) connected  
 * @remark the elements are dim-1 connected if for each pair of elements there
 *         exists a path via dim-1 adjacencies between them
 *
 * @param mesh (In) the mesh
 * @param rgns (In) rgns to check for connected-ness
 * @return true if connected, false o.w. 
 */
bool isFaceConnected(Mesh*& mesh, eArr& rgns) {
   assert( rgns.getSize() );
   set<MeshEntity*> visited;
   set<MeshEntity*> cavity;

   APF_ITERATE(eArr, rgns, rit) 
      cavity.insert(*rit);

   size_t count = 0;
   eList elms;
   MeshEntity* elm;
   // start the walk
   elms.push_back(rgns[0]); 
   while( ! elms.empty() ) {
      elm = elms.front();
      elms.pop_front();
      assert(elm);
      if( 1 == visited.count(elm) ) continue;
      visited.insert(elm);
      count++;
      eArr adjElms;
      getDwn2ndAdj(mesh, elm, adjElms);
      APF_ITERATE(eArr, adjElms, eit) {
         if ( 0 == visited.count(*eit) && 1 == cavity.count(*eit) )
            elms.push_back(*eit);
      }
      if ( count > rgns.getSize() ) {
        error("[%d] count > cavity size %zu > %ld\n", 
            mesh->getId(), count, (long)(rgns.getSize()));
        exit(EXIT_FAILURE);
      }
   }

   if ( count == rgns.getSize() ) {
      return true;
   } else {
      return false;
   }
}

inline int getCavityHighSurfAdjPid(Mesh* mesh, MeshTag* vtag, eArr& rgns) {
   set<MeshEntity*> edges;
   //eArr adjEdges;
   Downward adjEdges;
   APF_ITERATE(eArr, rgns, eit) {
     if ( mesh->hasTag(*eit, vtag) ) continue;
     const int na = mesh->getDownward(*eit, 1, adjEdges);
     for(int i=0; i<na; i++) 
       edges.insert(adjEdges[i]); 
   }
   if ( edges.size() == 0 )
     return -1;

   //  adjPartId, numShared edges
   typedef std::map<int, int> mapInt;
   mapInt sharedEdgeCount;
   APF_ITERATE(set<MeshEntity*>, edges, edgeItr) {
     apf::Copies rmt;
     mesh->getRemotes(*edgeItr, rmt);
     APF_ITERATE(apf::Copies, rmt, cpyItr) 
       ++(sharedEdgeCount[cpyItr->first]);
   }
   int maxPid = -1;
   int maxCnt = -1;
   APF_ITERATE(mapInt, sharedEdgeCount, cntItr) {
     if ( cntItr->second > maxCnt ) {
       maxPid = cntItr->first;
       maxCnt = cntItr->second;
     }
   }
   assert(-1 != maxPid);
   return maxPid;
}

inline bool isCavityTagged(Mesh* mesh, MeshTag* vtag, eArr& rgns) {
   APF_ITERATE(eArr, rgns, eit) 
      if ( mesh->hasTag(*eit, vtag) ) 
        return true;
   return false;
}

inline void printTotalMigr(const int itr, double w) {
   PCU_Add_Doubles(&w, 1);
   if( 0 == PCU_Comm_Self() ) {
     if( 0 == itr )
       debug(true, "migrInfo itr improveEntBalWeight\n");
     debug(true, "migrInfo %d %.3f\n", itr, w);
   }
}

inline bool isStuck(double weightMigr, const int entDim, 
    const int itr, const int dbgLvl) {
  double maxMigr = weightMigr;
  PCU_Max_Doubles(&maxMigr, 1);
  if ( isEql(maxMigr, 0) ) {
    if ( 0 == PCU_Comm_Self() )
      debug(dbgLvl>0,"No elements migrated in round %d of ent"
          "type %d improvement\n", itr, entDim);
    return true;
  }
  return false;
}

inline bool isDone(imbInfo& imb, const int itr, const int entDim, 
    const double maxImb) {
  if ( isLess(imb.maxImb[entDim],maxImb) ) {
    if ( 0 == PCU_Comm_Self() )
      status("MaxImb of ent type %d reached in round %d\n", entDim, itr);
    return true;
  }
  return false;
}

} //end unnamed namespace

double Parma::tagPtnMdlEdgeCavities(partInfo& part, const double maxW, 
    const size_t maxAdjElm, Migration* plan) {
  MeshEntity* vtx;
  eArr adjElms;
  eArr adjEdges;
  const int dim = mesh->getDimension();
  int planW = 0;
  int count = 0;

  MeshIterator* itr = mesh->begin(0);
  while( (vtx = mesh->iterate(itr)) ) {
    if( planW > maxW ) break;
    apf::Copies rmt;
    mesh->getRemotes(vtx, rmt);
    //  <><><><><> experimental code for migr across ptnMdl edges <><><><><><><>
    if( 1 < rmt.size() && mesh->isOwned(vtx) ) {
      mesh->getAdjacent(vtx, dim, adjElms);
      if( adjElms.getSize() < maxAdjElm  
          && ! isCavityTagged(mesh, vtag, adjElms)
          && isFaceConnected(mesh, adjElms) ) {
        int destPid = getCavityHighSurfAdjPid(mesh, vtag, adjElms);
        if( part.isCandidate(destPid) ) {
          planW++;
          count++;
          APF_ITERATE(eArr, adjElms, eit) {
            assert ( ! mesh->hasTag(*eit, vtag) );
            mesh->setIntTag(*eit, vtag, &destPid); 
            plan->send(*eit, destPid);
          }
        }
      }
    } 
  }
  mesh->end(itr);

  if ( planW > 0 ) 
    debug(dbgLvl>1, "[%d] %s tagged pMdlEdge %d\n", PCU_Comm_Self(), __func__, count);

  return planW;
}


/**
 * @brief compute the migration plan by finding vertices on the part boundary
 *        bounding bounding a small number of elements
 *
 * @param part (In) source part to migrate elements from
 * @param maxW (In) maximum weight to migrate
 * @param plan (InOut) migration plan for mesh elements 
 * @return sum of element weights in the plan
 */
double Parma::tagElmsForMigr_Vtx(partInfo& part, const double maxW, 
    Migration* plan) {
   const int maxBoundedElm = 6;
   double planW=0;
   for( int maxAdjElm=2; maxAdjElm<=maxBoundedElm; maxAdjElm+=2) 
     planW += tagSmallCavitiesForMigr(part, maxW-planW, maxAdjElm, plan);

   for( int maxAdjElm=2; maxAdjElm<=maxBoundedElm; maxAdjElm+=2) 
     planW += tagPtnMdlEdgeCavities(part, maxW-planW, maxAdjElm, plan);

   return planW;
}

double Parma::tagElmsForMigr(partInfo& part, const double maxW, 
    Migration* plan) {
   const int maxFaceOnPb = 4;
   const int minFaceOnPb = 2;
   double planW=0;
   for( int maxface=maxFaceOnPb; maxface>=minFaceOnPb; maxface--) 
     planW += tagElmsForMigr(part, maxW-planW, maxface, plan);

   return planW;
}

/**
 * @brief compute the migration plan by finding elements with high surface 
 *        area on the part boundary 
 *
 * @param part (In) source part to migrate elements from
 * @param plan (InOut) migration plan for mesh elements 
 * @return sum of element weights in the plan
 */
double Parma::tagElmsForMigr(partInfo& part, const double maxW, 
    const int maxFace, Migration* plan){

   double planWeight = 0;

   MeshEntity* face;
   MeshEntity* upElm;
   Downward adjFaces;

   const int dim = mesh->getDimension();
   MeshIterator* itr = mesh->begin(dim-1);
   while( (face = mesh->iterate(itr)) ) {
     if( planWeight > maxW ) break;
     apf::Copies rmt;
     mesh->getRemotes(face, rmt);
     if ( 1 != rmt.size()) continue;
     int destPid = (rmt.begin())->first;
     if( part.isCandidate(destPid) ) {
       upElm = getUpEnt(mesh, face);
       const int na = mesh->getDownward(upElm, dim-1, adjFaces);
       const int numPbFaces = getNumFaceOnPb(mesh, destPid, adjFaces, na);
       if( numPbFaces >= maxFace && !mesh->hasTag(upElm, vtag) ) {
         planWeight += addToPlan(mesh, wtag, upElm, destPid, plan);
         markAsVisited(mesh, vtag, destPid, upElm, adjFaces, na);
       }
     }
   }
   mesh->end(itr);
   return planWeight;
}

double Parma::tagSmallCavitiesForMigr(partInfo& part, const double maxW, 
    const size_t maxAdjElm, Migration* plan) {
  int planWeight = 0;

  MeshEntity* vtx;
  eArr adjElms;
  const int dim = mesh->getDimension();

  MeshIterator* itr = mesh->begin(0);
  while( (vtx = mesh->iterate(itr)) ) {
    if( planWeight > maxW ) break;
    apf::Copies rmt;
    mesh->getRemotes(vtx, rmt);
    if( 1 == rmt.size() ) {
      int destPid = (rmt.begin())->first;
      if( part.isCandidate(destPid) ) {
        mesh->getAdjacent(vtx, dim, adjElms);
        if( adjElms.getSize() < maxAdjElm ) {
          planWeight++;
          APF_ITERATE(eArr, adjElms, eit) {
            if ( mesh->hasTag(*eit, vtag) ) continue;
            mesh->setIntTag(*eit, vtag, &destPid); 
            plan->send(*eit, destPid);
          }
        }
      }
    }
  }
  mesh->end(itr);

  return planWeight;
}


/**
 * @brief improve the entity balance for the specified entity type
 *
 * @param pl (In) priority list
 * @param plIdx (In) priority list idx
 * @param maxImb (In) maximum entity imbalance > 1.0
 * @param itr (In) iteration number
 * @return weight of elements migrated on this process
 */
double Parma::improveEntBalance(const priorityList &pl, const int plIdx, 
      Array<double,4> avgW, const int) {
   if ( pl.priority[plIdx] <= 0 ) return 0;

   const int entDim = pl.entDim[plIdx];

   double maxImbW[4];
   for (int i=0;i<4;i++) maxImbW[i] = avgW[i]*maxImb;
   partInfo part(mesh, wtag, pl, plIdx, maxImbW);
   if( entDim != 0 ) 
     part.getGreedyMigrationSchedule2(pl, plIdx, maxImbW);

   const double maxWeight = part.weight[entDim] - avgW[entDim];
   Migration* plan = new Migration(mesh);
   double totW = 0;
   if( entDim == 0 ) 
     totW = tagElmsForMigr_Vtx(part, maxWeight, plan);
   else 
     totW = tagElmsForMigr(part, maxWeight, plan);

   clearTag(mesh, vtag);

   for(int i=0; i<plan->count(); i++)
     assert( 3 == apf::getDimension(mesh, plan->get(i)) );
   mesh->migrate(plan); //plan deleted here
   return totW;
}


void Parma::improveBalance(const priorityList& pl, const int plIdx) {
  const int entDim = pl.entDim[plIdx];
  if ( pl.priority[plIdx] <= 0 ) return;
  if ( 0 == PCU_Comm_Self() ) 
    debug(dbgLvl>0, "Improving ent dimension %d balance\n", entDim); 

  imbInfo imb;
  int itr=0;
  for (; itr < maxIter; itr++) {
    imb.get(mesh,wtag);
    if( dbgLvl>0 ) imb.print();
    if( isDone(imb, itr, entDim, maxImb) ) break;

    double wMigr = improveEntBalance(pl, plIdx, imb.minAvgW, itr);
    imb.get(mesh,wtag);
    if( dbgLvl>0 ) printTotalMigr(itr, wMigr);
    if( isStuck(wMigr, entDim, itr, dbgLvl) ) break;
    if( itr == maxIter-1 && 0 == PCU_Comm_Self() ) 
      debug(dbgLvl>0,"Max iterations reached for ent dimension %d\n", entDim); 
  }
}


/** 
 * @brief run partition improvement 
 * @remark see http://redmine.scorec.rpi.edu/projects/parma/wiki for more info
 *
 * @param mesh (InOut) partitioned mesh
 * @param priority (In) four integers representing the priority of mesh entity types to be balanced [vtx,edge,face,rgn]
 * @param dbgLvl (In) 0: off, >0 increasing amounts of runtime information
 * @param maxIter (In) maximum number of improvement iterations per entity type
 * @param maxImb (In) maximum entity imbalance tolerance > 1.0
 *
 * @return zero on success, non-zero otherwise
 */ 
int Parma::run(int (*priority)[4], const int dbgLvl_, 
    const int maxIter_, const double maxImb_) { 
   //ParMA_Embed_Version();
   if ( 0 == PCU_Comm_Self() )
     ParMA_Print_Version();

   if( PCU_Comm_Peers() == 1 ) return 0; // can't do much with a serial mesh

   this->dbgLvl = dbgLvl_;
   this->maxIter = maxIter_;
   this->maxImb = maxImb_;

   double t1 = MPI_Wtime();
   priorityList pl;
   pl.sort(priority);
   if( 0 == PCU_Comm_Self() ) pl.print();

   for (int plIdx=0; plIdx < 4; plIdx++)
     improveBalance(pl, plIdx);
   printElapsedTime(__func__, MPI_Wtime() - t1);

   return 0;
}


Parma::Parma(Mesh*& m) {
   init(m);
   userWeights = false;
   wtag = mesh->createDoubleTag("weight",1);
   assert(wtag);
   apf::MeshEntity* e;
   const int dim = mesh->getDimension();
   apf::Downward adj;
   MeshIterator* itr = mesh->begin(dim);
   while( (e = mesh->iterate(itr)) ) {
     int nd = mesh->getDownward(e, 0, adj);
     assert(nd>0);
     double w = static_cast<double>(nd);
     mesh->setDoubleTag(e, wtag, &w);
   }
   mesh->end(itr);

   vtag = mesh->createIntTag("visited",1);
}

Parma::Parma(Mesh*& m, MeshTag* t) {
   init(m);
   userWeights = true;
   wtag = t;
   vtag = mesh->createIntTag("visited",1);
}

void Parma::init(Mesh*& m) {
  mesh = m;
  if ( 3 != mesh->getDimension() ) {
    if ( 0 == PCU_Comm_Self() ) 
      error("ParMA supports 3D meshes only! "
          "input mesh is %dD ... exiting\n", 
          mesh->getDimension());
    exit(EXIT_FAILURE);
  }
}

Parma::~Parma() {
   MeshEntity* e;
   for (int d=0; d<=mesh->getDimension();d++) {
      MeshIterator* itr = mesh->begin(d);
      while( (e = mesh->iterate(itr)) ) {
         if( !userWeights && mesh->hasTag(e, wtag) ) {
            mesh->removeTag(e, wtag);
         }
         if( mesh->hasTag(e, vtag) ) {
            mesh->removeTag(e, vtag);
         }
      }
      mesh->end(itr);
   }
   if( !userWeights )
     mesh->destroyTag(wtag); 
   mesh->destroyTag(vtag); 
}
