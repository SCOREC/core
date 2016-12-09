#include "sam.h"
#include "samSz.h"
#include <apfField.h>
#include <cassert>
#include <stdio.h>

namespace sam {

apf::Field* errorThreshold(apf::Mesh* m, const char* fieldName,
    const unsigned idx, const double limit, const double factor)
{
  apf::Field* f = m->findField(fieldName);
  assert(f);
  apf::synchronize(f);
  apf::Field* newSz = apf::createFieldOn(m,"specifiedIso",apf::SCALAR);
  apf::Field* curSz = samSz::isoSize(m);
  double* vals = new double[f->countComponents()];
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    apf::getComponents(f, vtx, 0, vals);
    double h = apf::getScalar(curSz,vtx,0);
    if ( vals[idx] > limit )
      h *= factor;
    apf::setScalar(newSz,vtx,0,h);
  }
  m->end(itr);
  apf::destroyField(curSz);
  delete [] vals;
  return newSz;
}

apf::Field* compareIsoSF(apf::Mesh* m, const char* desiredSzFld, int method)
{
  apf::Field* f = m->findField(desiredSzFld);
  assert(f);
  assert(f->countComponents() == 1);
  apf::synchronize(f);
  apf::Field* newSz = apf::createFieldOn(m,"compareIsoSF",apf::SCALAR);
  apf::Field* curSz = samSz::isoSize(m);
  double dsr = 0.0;
  double cur = 0.0;
  double rsl = 0.0;
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    dsr = apf::getScalar(f, vtx, 0);
    cur = apf::getScalar(curSz,vtx,0);
    switch (method) {
      case 1:  rsl = (cur-dsr) / dsr; break; // ratio of difference
      default: rsl = cur / dsr;
    }
    apf::setScalar(newSz,vtx,0,rsl);
  }
  m->end(itr);
  apf::destroyField(curSz);
  return newSz;
}

apf::Field* specifiedIso(apf::Mesh* m, const char* fieldName, const unsigned idx)
{
  apf::Field* f = m->findField(fieldName);
  assert(f);
  apf::synchronize(f);
  apf::Field* newSz = apf::createFieldOn(m,"specifiedIso",apf::SCALAR);
  double* vals = new double[f->countComponents()];
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    apf::getComponents(f, vtx, 0, vals);
    double h = vals[idx];
    apf::setScalar(newSz,vtx,0,h);
  }
  m->end(itr);
  delete [] vals;
  return newSz;
}

static bool withinBox(apf::Vector3 points, double* box) {
  if ( fabs(points[0]- box[0]) <= box[3] &&
       fabs(points[1]- box[1]) <= box[4] &&
       fabs(points[2]- box[2]) <= box[5] )
    return true;
  else
    return false;
}

apf::Field* multiplySF(apf::Mesh* m, apf::Field* sf, double factor) {
  apf::Field* newSz = createFieldOn(m, "multiplySz", apf::SCALAR);
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    double h = apf::getScalar(sf,vtx,0);
    apf::setScalar(newSz,vtx,0,h*factor);
  }
  m->end(itr);
  apf::destroyField(sf);
  return newSz;
}

apf::Field* multiplySFBox(apf::Mesh* m, apf::Field* sf, double factor, double* box) {
  assert(box[3] > 0 && box[4] > 0 && box[5] > 0);
  apf::Field* newSz = apf::createFieldOn(m,"multiplySzBox",apf::SCALAR);
  apf::Vector3 points;
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    double h = apf::getScalar(sf,vtx,0);
    m->getPoint(vtx, 0, points);
    if (withinBox(points, box))
      apf::setScalar(newSz,vtx,0,h*factor);
    else
      apf::setScalar(newSz,vtx,0,h);
  }
  m->end(itr);
  apf::destroyField(sf);
  return newSz;
}

}
