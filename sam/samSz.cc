#include "samSz.h"
#include <apf.h>
#include <apfField.h>
#include <pcu_util.h>

namespace {

void getEdgeLenAndCnt(apf::Mesh* m, apf::MeshEntity* v, double& len, int& cnt) {
  len=0; cnt=0;
  apf::Up edges;
  m->getUp(v, edges);
  for(int eIdx=0; eIdx < edges.n; eIdx++) {
    if( m->isOwned(edges.e[eIdx])) {
      cnt++;
      len += measure(m,edges.e[eIdx]);
    }
  }
}

void setLenAndCnt(apf::Mesh* m, apf::Field* fLen, apf::Field* fCnt) {
  double len = 0;
  int cnt = 0;
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    getEdgeLenAndCnt(m, vtx, len, cnt);
    apf::setScalar(fLen, vtx, 0, len);
    apf::setScalar(fCnt, vtx, 0, cnt);
  }
  m->end(itr);
  apf::accumulate(fLen);
  apf::accumulate(fCnt);
  apf::synchronize(fLen);
  apf::synchronize(fCnt);
}

apf::Field* getIsoSize(apf::Mesh* m, apf::Field* fLen, apf::Field* fCnt) {
  if (m->findField("isoSize"))
    apf::destroyField(m->findField("isoSize"));
  apf::Field* sz = createFieldOn(m, "isoSize", apf::SCALAR);
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    int cnt = static_cast<int>(apf::getScalar(fCnt,vtx,0));
    double h = apf::getScalar(fLen,vtx,0);
    apf::setScalar(sz,vtx,0,h/cnt);
  }
  m->end(itr);
  return sz;
}

}

namespace samSz {

apf::Field* isoSize(apf::Mesh* m) {
  apf::Field* fLen = createFieldOn(m, "incidentEdgeLength", apf::SCALAR);
  apf::Field* fCnt = createFieldOn(m, "incidentEdgeCount", apf::SCALAR);
  setLenAndCnt(m,fLen,fCnt);
  apf::Field* isoSz = getIsoSize(m,fLen,fCnt);
  apf::destroyField(fLen);
  apf::destroyField(fCnt);
  return isoSz;
}

}
