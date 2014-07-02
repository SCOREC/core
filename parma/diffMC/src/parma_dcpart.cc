#include "parma_dcpart.h"
#include "parma_commons.h"
#include "PCU.h"
#include <stdio.h>
#include <vector>
#include <set>
#include <list>
#include <map>
#include "mpi.h"

using apf::DynamicArray;
using apf::Array;
using apf::MeshEntity;
using apf::Mesh;
using apf::MeshTag;
using apf::MeshIterator;
using apf::Migration;

using std::set;
using std::map;
using std::vector;

using parmaCommons::printElapsedTime;
using parmaCommons::debug;
using parmaCommons::error;

typedef std::list<MeshEntity*> eList;
typedef DynamicArray<MeshEntity*> eArr;

dcPart::dcPart() {
   exit(EXIT_FAILURE);
}

dcPart::dcPart(Mesh*& mesh) : m(mesh) {
   init(m);
}

void dcPart::init(Mesh*& mesh) {
   vtag = mesh->createIntTag("dcVisited",1);
}

inline void clearTag(Mesh*& m, MeshTag* t) {
   MeshEntity* e;
   MeshIterator* itr = m->begin(m->getDimension());
   while( (e = m->iterate(itr)) ) {
      if( m->hasTag(e, t) ) 
         m->removeTag(e, t);
   }
   m->end(itr);
}

dcPart::~dcPart() {
   clearTag(m, vtag);
   m->destroyTag(vtag); 
}

inline int getEntDim(Mesh* m, MeshEntity* e) {
   const int t = m->getType(e);
   return apf::Mesh::getEntityDimension(t);
}

inline MeshEntity* getUpElm(Mesh* m, MeshEntity* e) {
   const int upDim = getEntDim(m, e) + 1;
   eArr adjEnt;
   m->getAdjacent(e, upDim, adjEnt);
   assert( NULL != adjEnt[0] );
   return adjEnt[0];
}

inline void getDwn2ndAdj(Mesh* m, MeshEntity* elm, eArr& adj) {
   const int dim = getEntDim(m, elm);
   eArr adjF;  
   m->getAdjacent(elm, dim-1, adjF);
   APF_ITERATE(eArr, adjF, fit) {
      eArr adjElms;  
      m->getAdjacent(*fit, dim, adjElms);
      APF_ITERATE(eArr, adjElms, eit) 
         if ( *eit != elm ) 
            adj.append(*eit);
   }
}

int dcPart::numDisconnectedComps() {
   double t1 = MPI_Wtime();
   dcCompSz.clear();
   clearTag(m, vtag);
   int numDc = 0;
   size_t count = 0;
   while( count != m->count(m->getDimension()) ) {
      dcCompSz.push_back( walkPart(numDc) );
      count += dcCompSz[numDc];
      numDc++;
   }
   printElapsedTime(__func__, MPI_Wtime() - t1);
   return numDc-1;
}

int dcPart::walkPart(int visited) {
   size_t count = 0;
   const int dim = m->getDimension();

   eList elms;
   MeshEntity* elm;
   // find an untagged element
   MeshIterator* itr = m->begin(dim);
   while( (elm = m->iterate(itr)) && m->hasTag(elm, vtag) );
   m->end(itr);
   // start the walk
   assert(elm);
   elms.push_back(elm); 
   while( ! elms.empty() ) {
      elm = elms.front();
      elms.pop_front();
      assert( elm != NULL );
      if ( m->hasTag(elm, vtag) ) continue;
      m->setIntTag(elm, vtag, &visited); 
      count++;
      eArr adjElms;
      getDwn2ndAdj(m, elm, adjElms);
      APF_ITERATE(eArr, adjElms, eit) {
         if ( ! m->hasTag(*eit, vtag) ) 
            elms.push_back(*eit);
      }
      if ( count > m->count(dim) ) { 
	 error("[%d] count > part size %d > %lu\n", 
	       m->getId(), count, m->count(dim));
         exit(EXIT_FAILURE);
      }
   }
   return count;
}

/**
 * @brief remove the disconnected set(s) of elements from the part
 * @remark migrate the disconnected set(s) of elements into the adjacent part that
 *         shares the most faces with the disconnected set of elements
 *         requires that the sets of elements forming disconnected components are 
 *         tagged
 */
void dcPart::fix() {
   double t1 = MPI_Wtime();

   // < dcComId , mergeTgtPid >
   migrTgt dcCompTgts; 

   int maxSz = -1;
   APF_ITERATE(vector<int>, dcCompSz, dc) 
      if( *dc > maxSz ) 
         maxSz = *dc; 
   
   for(size_t i=0; i<dcCompSz.size(); i++)
      if( dcCompSz[i] != maxSz ) 
         dcCompTgts[i] = checkResidence(i);
   assert( dcCompTgts.size() == dcCompSz.size()-1 );
   Migration* plan = new Migration(m);
   setupPlan(dcCompTgts, plan);
   clearTag(m, vtag);
   m->migrate(plan); // plan deleted here
   printElapsedTime(__func__, MPI_Wtime() - t1);
}

int dcPart::checkResidence(const int dcComp) {
   // < dcComId, maxFace >
   typedef map<int, int> mii; 
   mii bdryFaceCnt; 

   const int dim = m->getDimension();
   int tval;
   eArr adjFaces;
   apf::Parts resPid;

   MeshEntity* e;
   MeshIterator* itr = m->begin(dim);
   while( (e = m->iterate(itr)) ) {
      if( m->hasTag(e, vtag) ) {
	 m->getIntTag(e, vtag, &tval);
	 if( tval == dcComp ) { 
	    m->getAdjacent(e, dim-1, adjFaces); 
	    APF_ITERATE(eArr, adjFaces, fit) {
               if( ! m->isShared(*fit) ) continue;
	       m->getResidence(*fit, resPid);
	       APF_ITERATE(apf::Parts, resPid, rp) {
                  if( bdryFaceCnt.count(*rp) ) {
		     bdryFaceCnt[*rp] += 1;
                  } else {
		     bdryFaceCnt[*rp] = 0;
                  }
	       }
	    }
	 }
      }
   }
   m->end(itr);
   int max = -1;
   int maxId = -1;
   APF_ITERATE(mii , bdryFaceCnt, bf) {
      if( bf->first != m->getId() && bf->second > max ) {
         max = bf->second;
         maxId = bf->first;
      }
      debug(false, "[%d] bdryFaceCnt dcComp %d ap %d faces %d\n", 
	    PCU_Comm_Self(), dcComp, bf->first, bf->second);
   }
   debug(false, "[%d] %s maxId %d max %d\n", 
	 PCU_Comm_Self(), __func__, maxId, max);
   return maxId;
}

void dcPart::setupPlan(migrTgt& dcCompTgts, Migration* plan) {
   MeshEntity* e;
   MeshIterator* itr = m->begin(m->getDimension());
   while( (e = m->iterate(itr)) ) {
      if( m->hasTag(e, vtag) ) {
	 int tval = -1;
	 m->getIntTag(e, vtag, &tval);
	 if ( dcCompTgts.count(tval) ) {
            plan->send(e, dcCompTgts[tval]);
	 }
      }
   }
   m->end(itr);
   debug(false, "[%d] fixDcComps migrating %d\n", m->getId(), plan->count());
}
 

inline bool isShared(Mesh* m, MeshEntity* elm) {
   const int dim = getEntDim(m, elm);
   assert( dim == m->getDimension() );
   eArr adjEnt;
   m->getAdjacent(elm, dim-1, adjEnt);
   APF_ITERATE(eArr, adjEnt, it) {
      if( m->isShared(*it) ) 
         return true;
   }
   return false;
}

void dcPart::makeDisconnectedComps(const int numDcComps) {
   clearTag(m, vtag);

   const int destPid = (m->getId() + 1) % PCU_Comm_Peers();

   Migration* plan = new Migration(m);
   for(int i=0; i<numDcComps; i++) {
      eList elms;
      MeshEntity* elm;
      bool found = false;
      // find an untagged element
      MeshIterator* itr = m->begin(m->getDimension());
      while( (elm = m->iterate(itr)) && !found ) {
         if ( m->hasTag(elm, vtag) ) continue;
	 eArr adjElms;
         // check if face adj elms are tagged
	 getDwn2ndAdj(m, elm, adjElms);
         int numDirtyElms = 0;
	 APF_ITERATE(eArr, adjElms, eit) 
	    if( m->hasTag(*eit, vtag) || isShared(m, *eit) )  
               numDirtyElms++;
         if ( numDirtyElms == 0 ) {
            plan->send(elm, destPid);
	    m->setIntTag(elm, vtag, &i); 
            found = true;
	 }
      }
      m->end(itr);
   }
   
   debug(false, "[%d] migrating %d to %d\n", 
	 m->getId(), plan->count(), destPid);
   m->migrate(plan); //plan deleted here
}
