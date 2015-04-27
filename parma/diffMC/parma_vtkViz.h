#ifndef PARMA_VTKVIZ_H
#define PARMA_VTKVIZ_H

#include <sstream>
#include <string>
#include <parma.h>
#include <PCU.h>
#include <apfShape.h>
#include <apfNumbering.h>
#include "parma_dcpart.h"

#define TOSIZET(a) static_cast<size_t>(a)
#define TOINT(a) static_cast<int>(a)
#define TODBL(a) static_cast<double>(a)

namespace parma {
  apf::Numbering* makeVtxNumbering(apf::Mesh* m, apf::MeshTag* t) {
    apf::Numbering* n = 
      apf::createNumbering(m,"parmaVtxNumbering",m->getShape(),1);
    apf::MeshEntity* v;
    apf::MeshIterator* it = m->begin(0);
    const int node = 0, comp = 0;
    int dist;
    while( (v = m->iterate(it)) ) {
      if( m->hasTag(v,t) )
        m->getIntTag(v,t,&dist);
      else 
        dist = 999999;
      apf::number(n,v,node,comp,dist);
    }
    m->end(it);
    return n;
  }

  void writeAllVtk(apf::Mesh* m, const char* key="parma") {
    static int stepCnt = 0;
    std::stringstream ss;
    ss << key << '.' << stepCnt << '.';
    std::string pre = ss.str();
    apf::writeVtkFiles(pre.c_str(), m);
    stepCnt++;
  }

  void writeVtk(apf::Mesh* m, const char* key, int step) {
    std::stringstream ss;
    ss << key << '.' << step << '.';
    std::string pre = ss.str();
    apf::writeOneVtkFile(pre.c_str(), m);
  }

  void writeMaxParts(apf::Mesh* m, const char* key="", parma::dcComponents* c=NULL) {
    static int stepCnt = 0;
    long tot;
    int min, max, loc; 
    double avg;
    if( c ) {
      int totDc;
      totDc = max = loc = TOINT(c->size());
      PCU_Max_Ints(&max, 1);
      PCU_Add_Ints(&totDc, 1);
      avg = TODBL(totDc)/PCU_Comm_Peers();
    } else {
      Parma_GetDisconnectedStats(m, max, avg, loc);
    }
    if( loc == max ) {
      std::string pre = std::string("maxDc") + std::string(key);
      writeVtk(m, pre.c_str(), stepCnt);
    }
    Parma_GetEntStats(m,0,tot,min,max,avg,loc);
    if( loc == max ) {
      std::string pre = std::string("maxVtxImb") + std::string(key);
      writeVtk(m, pre.c_str(), stepCnt);
    }
    stepCnt++;
  }
}
#endif
