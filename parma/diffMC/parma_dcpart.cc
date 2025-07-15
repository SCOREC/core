#include "parma_dcpart.h"
#include "parma_commons.h"
#include "parma_meshaux.h"
#include "parma_convert.h"
#include <maximalIndependentSet/mis.h>
#include <stdio.h>
#include <set>
#include <list>
#include <map>
#include <pcu_util.h>

typedef std::map<unsigned, unsigned> muu;

dcPart::dcPart(apf::Mesh*& mesh, unsigned v)
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

apf::MeshEntity* dcPart::getSeedEnt(unsigned i) {
  apf::MeshEntity* elm;
  apf::MeshIterator* itr = m->begin(m->getDimension());
  while( (elm = m->iterate(itr)) )
    if( ! isIsolated(elm) && (compId(elm) == i ) )
      break;
  m->end(itr);
  return elm;
}

unsigned dcPart::compId(apf::MeshEntity* e) {
  PCU_ALWAYS_ASSERT(m->hasTag(e, vtag));
  int c; m->getIntTag(e, vtag, &c);
  return TO_UINT(c);
}

bool dcPart::isIsolated(apf::MeshEntity* e) {
  return m->hasTag(e, isotag);
}

void dcPart::markIsolated(const unsigned dcComp) {
  int one = 1;
  int tval = -1;
  apf::MeshEntity* elm;
  apf::MeshIterator* itr = m->begin(m->getDimension());
  while( (elm = m->iterate(itr)) ) {
    if( m->hasTag(elm, vtag) ) {
      m->getIntTag(elm, vtag, &tval);
      if( tval == TO_INT(dcComp) ) {
        m->removeTag(elm, vtag); //clear the dc comp id
        m->setIntTag(elm, isotag, &one);
      }
    }
  }
  m->end(itr);
}

unsigned dcPart::getNumComps() {
  return TO_UINT(dcCompSz.size());
}

unsigned dcPart::getNumDcComps() {
  unsigned c = TO_UINT(dcCompSz.size());
  //If there are c components, and c > 0, then the 0th component
  //is considered the core so subtract one from c to get the number
  //of disconnected components.
  if( c > 0 ) c-=1;
  return c;
}

unsigned dcPart::getNumIso() {
  return numIso;
}

unsigned dcPart::getCompSize(unsigned i) {
  PCU_ALWAYS_ASSERT( i < dcCompSz.size());
  return dcCompSz[i];
}

unsigned dcPart::getCompPeer(unsigned i) {
  PCU_ALWAYS_ASSERT( i < dcCompSz.size());
  return dcCompNbor[i];
}

void dcPart::reset() {
   dcCompSz.clear();
   dcCompNbor.clear();
   clearTag(m, vtag);
   clearTag(m, isotag);
   numIso = 0;
}

unsigned dcPart::numDisconnectedComps() {
   double t1 = pcu::Time();
   reset();
   unsigned numDc = 0;
   size_t count = 0;
   unsigned self = TO_UINT(m->getId());
   const size_t numElms = m->count(m->getDimension());
   while( count != numElms ) {
      unsigned sz = walkPart(numDc);
      unsigned nbor = maxContactNeighbor(numDc);
      if( nbor != self || m->getPCU()->Peers() == 1 ) {
        dcCompSz.push_back(sz);
        dcCompNbor.push_back(nbor);
        numDc++;
      } else {
        numIso++;
        markIsolated(numDc);
      }
      count += sz;
   }
   if( verbose )
     parmaCommons::printElapsedTime(__func__, pcu::Time() - t1, m->getPCU());
   PCU_ALWAYS_ASSERT(numDc+numIso >= 1);
   return (numDc+numIso)-1;
}

unsigned dcPart::walkPart(unsigned visited) {
   unsigned count = 0;
   const int dim = m->getDimension();

   std::list<apf::MeshEntity*> elms;
   apf::MeshEntity* elm;
   // find an untagged element
   apf::MeshIterator* itr = m->begin(dim);
   while( (elm = m->iterate(itr)) &&
          ( m->hasTag(elm, vtag) || m->hasTag(elm, isotag) ) );
   m->end(itr);
   // start the walk
   PCU_ALWAYS_ASSERT(elm);
   elms.push_back(elm);
   while( ! elms.empty() ) {
      elm = elms.front();
      elms.pop_front();
      PCU_ALWAYS_ASSERT( elm != NULL );
      if ( m->hasTag(elm, vtag) || m->hasTag(elm, isotag) ) continue;
      const int dcId = TO_INT(visited);
      m->setIntTag(elm, vtag, &dcId);
      count++;
      apf::Adjacent adjElms;
      getDwn2ndAdj(m, elm, adjElms);
      APF_ITERATE(apf::Adjacent, adjElms, eit)
        if ( ! m->hasTag(*eit, vtag) )
          elms.push_back(*eit);
      PCU_ALWAYS_ASSERT( count <= TO_UINT(m->count(dim)) );
   }
   return count;
}


unsigned dcPart::maxContactNeighbor(const unsigned dcComp) {
   // < dcComId, maxFace >
   muu bdryFaceCnt;

   const int dim = m->getDimension();
   int tval;
   apf::Downward sides;
   apf::Parts resPid;

   apf::MeshEntity* e;
   apf::MeshIterator* itr = m->begin(dim);
   while( (e = m->iterate(itr)) ) {
      if( m->hasTag(e, vtag) ) {
	 m->getIntTag(e, vtag, &tval);
	 if( tval != TO_INT(dcComp) ) continue;
         int ns = m->getDownward(e, dim-1, sides);
         for(int sIdx=0; sIdx<ns; sIdx++) {
           apf::MeshEntity* s = sides[sIdx];
           if( ! m->isShared(s) ) continue;
           m->getResidence(s, resPid);
           APF_ITERATE(apf::Parts, resPid, rp)
             (bdryFaceCnt[TO_UINT(*rp)])++;
         }
      }
   }
   m->end(itr);
   unsigned max = 0;
   unsigned maxId = TO_UINT(m->getId());
   unsigned self = maxId;
   APF_ITERATE(muu, bdryFaceCnt, bf) {
      if( bf->first != self && bf->second > max ) {
         max = bf->second;
         maxId = bf->first;
      }
   }
   return maxId;
}
