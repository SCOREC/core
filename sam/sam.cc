#include "sam.h"
#include "samSz.h"
#include <apfField.h>
#include <cassert>
#include <stdio.h>
#include <math.h>

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

static bool withinCyl(apf::Vector3 points, double* cyl) {
  double length = sqrt(cyl[3]*cyl[3]+cyl[4]*cyl[4]+cyl[5]*cyl[5]);
  cyl[3] = cyl[3] / length;
  cyl[4] = cyl[4] / length;
  cyl[5] = cyl[5] / length;
  if((fabs((points[0]- cyl[0])*cyl[3]+
           (points[1]- cyl[1])*cyl[4]+
		   (points[2]- cyl[2])*cyl[5]) <= cyl[6]) &&
  sqrt(((points[0]- cyl[0])*cyl[4] - (points[1]- cyl[1])*cyl[3])*
       ((points[0]- cyl[0])*cyl[4] - (points[1]- cyl[1])*cyl[3])+
	   ((points[0]- cyl[0])*cyl[5] - (points[2]- cyl[2])*cyl[3])*
	   ((points[0]- cyl[0])*cyl[5] - (points[2]- cyl[2])*cyl[3])+
	   ((points[1]- cyl[1])*cyl[5] - (points[2]- cyl[2])*cyl[4])*
	   ((points[1]- cyl[1])*cyl[5] - (points[2]- cyl[2])*cyl[4]) <= cyl[7]))
    return true;
  else
    return false;
}

void multiplySF(apf::Mesh* m, apf::Field* sf, double factor) {
  assert(m->findField(apf::getName(sf)));
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    double h = apf::getScalar(sf,vtx,0);
    apf::setScalar(sf,vtx,0,h*factor);
  }
  m->end(itr);
}

void multiplySFBox(apf::Mesh* m, apf::Field* sf, double factor, double* box) {
  assert(box[3] > 0 && box[4] > 0 && box[5] > 0);
  assert(m->findField(apf::getName(sf)));
  apf::Vector3 points;
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    m->getPoint(vtx, 0, points);
    if (withinBox(points, box)) {
      double h = apf::getScalar(sf,vtx,0);
      apf::setScalar(sf,vtx,0,h*factor);
	}
  }
  m->end(itr);
}

void multiplySFCyl(apf::Mesh* m, apf::Field* sf, double factor, double* cyl) {
  /* cylinder is defined as {center_x, center_y, center_z,
     normal_x, normal_y, normal_z, half_height, radius}   */
  assert(cyl[6] > 0 && cyl[7] > 0);
  assert(m->findField(apf::getName(sf)));
  apf::Vector3 points;
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    m->getPoint(vtx, 0, points);
    if (withinCyl(points, cyl)) {
      double h = apf::getScalar(sf,vtx,0);
      apf::setScalar(sf,vtx,0,h*factor);
	}
  }
  m->end(itr);
}

void multiplySFRegion(apf::Mesh* m, apf::Field* sf, double factor, int tag) {
  /* tag should be a tag of a region */
  assert(m->findField(apf::getName(sf)));
  int type = 3; // 3D region by default
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    apf::ModelEntity* m_vtx = m->toModel(vtx);
    apf::ModelEntity* region = m->findModelEntity(type, tag);
    if (m_vtx == region ||
	    m->isInClosureOf(m_vtx, region)) {
      double h = apf::getScalar(sf,vtx,0);
      apf::setScalar(sf,vtx,0,h*factor);
	}
  }
  m->end(itr);
}

}
