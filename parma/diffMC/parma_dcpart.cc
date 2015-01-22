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
typedef map<unsigned, unsigned> muu;

namespace {
  bool isInMis(muu& mt) {
    unsigned seed = PCU_Comm_Self()+1;
    mis_init(seed);
    misLuby::partInfo part;
    part.id = PCU_Comm_Self();
    std::set<int> targets;
    APF_ITERATE(muu, mt, mtItr) {
      if( !targets.count(mtItr->second) ) {
        part.adjPartIds.push_back(mtItr->second);
        part.net.push_back(mtItr->second);
        targets.insert(mtItr->second);
      }
    }
    part.net.push_back(part.id);
    return mis(part,false,true);
  }

  inline void getDwn2ndAdj(apf::Mesh* m, apf::MeshEntity* ent, eArr& adj) {
    const int dim = apf::getDimension(m, ent);
    apf::getBridgeAdjacent(m, ent, dim - 1, dim, adj);
  }
}

dcPart::dcPart(Mesh*& mesh, unsigned v)
  : m(mesh), verbose(v) {
   vtag = mesh->createIntTag("dcVisited",1);
   isotag = mesh->createIntTag("dcIsolated",1);
   numDisconnectedComps();
}

dcPart::~dcPart() {
   clearTag(m, vtag);
   m->destroyTag(vtag);

   clearTag(m, isotag);
   m->destroyTag(isotag);
}

bool dcPart::isIsolated(apf::MeshEntity* e) {
  return m->hasTag(e, isotag);
}

void dcPart::markIsolated(const unsigned dcComp) {
  int one = 1;
  int tval = -1;
  MeshEntity* elm;
  MeshIterator* itr = m->begin(m->getDimension());
  while( (elm = m->iterate(itr)) ) {
    if( m->hasTag(elm, vtag) ) {
      m->getIntTag(elm, vtag, &tval);
      if( tval == static_cast<int>(dcComp) )
        m->setIntTag(elm, isotag, &one);
    }
  }
  m->end(itr);
}

unsigned dcPart::numDisconnectedComps() {
   double t1 = PCU_Time();
   dcCompSz.clear();
   dcCompNbor.clear();
   clearTag(m, vtag);
   clearTag(m, isotag);
   unsigned numDc = 0;
   unsigned count = 0;
   unsigned self = m->getId();
   const int dim = m->getDimension();
   const unsigned numElms = m->count(dim);
   while( count != numElms ) {
      unsigned sz = walkPart(numDc);
      unsigned nbor = maxContactNeighbor(numDc);
      if( nbor != self ) {
        dcCompSz.push_back(sz);
        dcCompNbor.push_back(nbor);
        numDc++;
      } else {
        markIsolated(numDc);
      }
      count += sz;
   }
   if( verbose )
     printElapsedTime(__func__, PCU_Time() - t1);
   assert(numDc >= 1);
   return numDc-1;
}

int dcPart::totNumDc() {
  int ndc = numDisconnectedComps();
  PCU_Add_Ints(&ndc, 1);
  return ndc;
}

unsigned dcPart::walkPart(unsigned visited) {
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
      const int dcId = visited;
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
    muu dcCompTgts;

    unsigned maxSz = 0;
    APF_ITERATE(vector<unsigned>, dcCompSz, dc)
      if( *dc > maxSz )
        maxSz = *dc;

    for(size_t i=0; i<dcCompSz.size(); i++)
      if( dcCompSz[i] != maxSz )
        dcCompTgts[i] = dcCompNbor[i];
    assert( dcCompTgts.size() == dcCompSz.size()-1 );
    Migration* plan = new Migration(m);
    if ( isInMis(dcCompTgts) )
      setupPlan(dcCompTgts, plan);

    clearTag(m, vtag);
    double t3 = PCU_Time();
    m->migrate(plan);
    if( 0 == PCU_Comm_Self() && verbose)
      status("loop %d components %d seconds <fix migrate> %.3f %.3f\n",
          loop, ndc, t3-t2, PCU_Time()-t3);
  }
  printElapsedTime(__func__, PCU_Time() - t1);
}

unsigned dcPart::maxContactNeighbor(const unsigned dcComp) {
   // < dcComId, maxFace >
   muu bdryFaceCnt;

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
   unsigned max = 0;
   unsigned maxId = m->getId();
   unsigned self = m->getId();
   APF_ITERATE(muu , bdryFaceCnt, bf) {
      if( bf->first != self && bf->second > max ) {
         max = bf->second;
         maxId = bf->first;
      }
   }
   return maxId;
}

void dcPart::setupPlan(muu& dcCompTgts, Migration* plan) {
   MeshEntity* e;
   MeshIterator* itr = m->begin(m->getDimension());
   while( (e = m->iterate(itr)) ) {
      if( m->hasTag(e, vtag) ) {
	 int tval = -1;
	 m->getIntTag(e, vtag, &tval);
         const unsigned dcId = tval;
	 if ( dcCompTgts.count(dcId) )
            plan->send(e, dcCompTgts[dcId]);
      }
   }
   m->end(itr);
}



