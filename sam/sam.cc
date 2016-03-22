#include "sam.h"
#include "samSz.h"
#include <cassert>
#include <stdio.h>

namespace sam {
  apf::Field* specifiedIso(apf::Mesh* m, const char* fieldName, 
      const unsigned idx, const double limit, const double factor) {
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
}
