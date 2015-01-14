#include "PCU.h"
#include "parma_dcpart.h"
#include "parma_commons.h"
#include "parma_meshaux.h"
#include <maximalIndependentSet/mis.h>
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
using parmaCommons::status;
using parmaCommons::error;

typedef std::list<MeshEntity*> eList;
typedef DynamicArray<MeshEntity*> eArr;

namespace {
  bool isInMis(migrTgt& mt) {
    unsigned int seed = static_cast<unsigned int>(PCU_Comm_Self()+1);
    mis_init(seed);
    misLuby::partInfo part;
    part.id = PCU_Comm_Self();
    std::set<int> targets;
    APF_ITERATE(migrTgt, mt, mtItr) {
      if( !targets.count(mtItr->second) ) {
        part.adjPartIds.push_back(mtItr->second);
        part.net.push_back(mtItr->second);
        targets.insert(mtItr->second);
      }
    }
    part.net.push_back(part.id);
    return mis(part,false,true);
  }
}

dcPart::dcPart(Mesh*& mesh) : m(mesh) {
   init(m);
}

void dcPart::init(Mesh*& mesh) {
   vtag = mesh->createIntTag("dcVisited",1);
}

dcPart::~dcPart() {
   clearTag(m, vtag);
   m->destroyTag(vtag);
}

inline MeshEntity* getUpElm(Mesh* m, MeshEntity* e) {
   const int upDim = apf::getDimension(m, e) + 1;
   eArr adjEnt;
   m->getAdjacent(e, upDim, adjEnt);
   assert( NULL != adjEnt[0] );
   return adjEnt[0];
}

int dcPart::numDisconnectedComps() {
   double t1 = PCU_Time();
   dcCompSz.clear();
   clearTag(m, vtag);
   size_t numDc = 0;
   int count = 0;
   const int dim = m->getDimension();
   const int numElms = static_cast<int>(m->count(dim));
   while( count != numElms ) {
      dcCompSz.push_back( walkPart(numDc) );
      count += dcCompSz[numDc];
      numDc++;
   }
   printElapsedTime(__func__, PCU_Time() - t1);
   return static_cast<int>(numDc-1);
}

int dcPart::totNumDc() {
  int ndc = numDisconnectedComps();
  PCU_Add_Ints(&ndc, 1);
  return ndc;
}

size_t dcPart::walkPart(size_t visited) {
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
      int dcId = static_cast<int>(visited);
      m->setIntTag(elm, vtag, &dcId);
      count++;
      eArr adjElms;
      getDwn2ndAdj(m, elm, adjElms);
      if (adjElms.getSize())
        APF_ITERATE(eArr, adjElms, eit) {
           if ( ! m->hasTag(*eit, vtag) )
              elms.push_back(*eit);
        }
      if ( count > m->count(dim) ) {
	 error("[%d] count > part size %d > %ld\n",
	       m->getId(), count, static_cast<long>(m->count(dim)));
         exit(EXIT_FAILURE);
      }
   }
   return count;
}

/**
 * @brief remove the disconnected set(s) of elements from the part
 * @remark migrate the disconnected set(s) of elements into the adjacent part
 *         that shares the most faces with the disconnected set of elements
 *         requires that the sets of elements forming disconnected components
 *         are tagged
 */
void dcPart::fix() {
  double t1 = PCU_Time();
  int loop = 0;
  int ndc = 0;
  while( (ndc = totNumDc()) && loop++ < 50 ) {
    double t2 = PCU_Time();
    migrTgt dcCompTgts;

    size_t maxSz = 0;
    APF_ITERATE(vector<size_t>, dcCompSz, dc)
      if( *dc > maxSz )
        maxSz = *dc;

    size_t isolated = 0;
    for(size_t i=0; i<dcCompSz.size(); i++)
      if( dcCompSz[i] != maxSz ) {
        int res = checkResidence(i);
        if ( res != -1 )
          dcCompTgts[i] = res;
        else
          isolated++;
      }
    assert( dcCompTgts.size() + isolated == dcCompSz.size()-1 );
    Migration* plan = new Migration(m);
    if ( isInMis(dcCompTgts) )
      setupPlan(dcCompTgts, plan);

    clearTag(m, vtag);
    double t3 = PCU_Time();
    m->migrate(plan);
    if( 0 == PCU_Comm_Self() )
      status("loop %d components %d seconds <fix migrate> %.3f %.3f\n",
          loop, ndc, t3-t2, PCU_Time()-t3);
  }
  printElapsedTime(__func__, PCU_Time() - t1);
}

int dcPart::checkResidence(const size_t dcComp) {
   // < dcComId, maxFace >
   typedef map<int, int> mii;
   mii bdryFaceCnt;

   const int dim = m->getDimension();
   int tval;
   apf::Downward sides;
   apf::Parts resPid;

   MeshEntity* e;
   MeshIterator* itr = m->begin(dim);
   while( (e = m->iterate(itr)) ) {
      if( m->hasTag(e, vtag) ) {
	 m->getIntTag(e, vtag, &tval);
	 if( tval != static_cast<int>(dcComp) ) continue;
         int ns = m->getDownward(e, dim-1, sides);
         for(int sIdx=0; sIdx<ns; sIdx++) {
           MeshEntity* s = sides[sIdx];
           if( ! m->isShared(s) ) continue;
           m->getResidence(s, resPid);
           APF_ITERATE(apf::Parts, resPid, rp)
             (bdryFaceCnt[*rp])++;
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
   }
   return maxId;
}

void dcPart::setupPlan(migrTgt& dcCompTgts, Migration* plan) {
   MeshEntity* e;
   MeshIterator* itr = m->begin(m->getDimension());
   while( (e = m->iterate(itr)) ) {
      if( m->hasTag(e, vtag) ) {
	 int tval = -1;
	 m->getIntTag(e, vtag, &tval);
         const size_t dcId = static_cast<size_t>(tval);
	 if ( dcCompTgts.count(dcId) ) {
            plan->send(e, dcCompTgts[dcId]);
	 }
      }
   }
   m->end(itr);
}


inline bool isShared(Mesh* m, MeshEntity* elm) {
   const int dim = apf::getDimension(m, elm);
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

   clearTag(m, vtag);
   m->migrate(plan); //plan deleted here
}
