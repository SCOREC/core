#include "sam.h"
#include "samSz.h"
#include <apfGeometry.h>
#include <apfField.h>
#include <pcu_util.h>
#include <stdio.h>
#include <math.h>

namespace sam {

apf::Field* errorThreshold(apf::Mesh* m, const char* fieldName,
    const unsigned idx, const double limit, const double factor)
{
  apf::Field* f = m->findField(fieldName);
  PCU_ALWAYS_ASSERT(f);
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
  PCU_ALWAYS_ASSERT(f);
  PCU_ALWAYS_ASSERT(f->countComponents() == 1);
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
  PCU_ALWAYS_ASSERT(f);
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

void multiplySF(apf::Mesh* m, apf::Field* sf, double factor) {
  PCU_ALWAYS_ASSERT(m->findField(apf::getName(sf)));
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    double h = apf::getScalar(sf,vtx,0);
    apf::setScalar(sf,vtx,0,h*factor);
  }
  m->end(itr);
}

void multiplySFBox(apf::Mesh* m, apf::Field* sf, double factor, double* box) {
  PCU_ALWAYS_ASSERT(box[3] > 0 && box[4] > 0 && box[5] > 0);
  PCU_ALWAYS_ASSERT(m->findField(apf::getName(sf)));
  apf::Vector3 points;
  apf::Vector3 center = apf::Vector3(box[0],box[1],box[2]);
  apf::Vector3 size   = apf::Vector3(box[3],box[4],box[5]);
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    m->getPoint(vtx, 0, points);
    if (apf::withinBox(points, center, size)) {
      double h = apf::getScalar(sf,vtx,0);
      apf::setScalar(sf,vtx,0,h*factor);
	}
  }
  m->end(itr);
}

void multiplySFCyl(apf::Mesh* m, apf::Field* sf, double factor, double* cyl) {
  /* cylinder is defined as {center_x, center_y, center_z,
     normal_x, normal_y, normal_z, half_height, radius}   */
  PCU_ALWAYS_ASSERT(cyl[6] > 0 && cyl[7] > 0);
  PCU_ALWAYS_ASSERT(m->findField(apf::getName(sf)));
  apf::Vector3 points;
  apf::Vector3 center = apf::Vector3(cyl[0],cyl[1],cyl[2]);
  apf::Vector3 normal = apf::Vector3(cyl[3],cyl[4],cyl[5]);
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    m->getPoint(vtx, 0, points);
    if (apf::withinCyl(points, center, cyl[6], cyl[7], normal)) {
      double h = apf::getScalar(sf,vtx,0);
      apf::setScalar(sf,vtx,0,h*factor);
	}
  }
  m->end(itr);
}

void multiplySFRegion(apf::Mesh* m, apf::Field* sf, double factor, int tag) {
  /* tag should be a tag of a region */
  PCU_ALWAYS_ASSERT(m->findField(apf::getName(sf)));
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
