#include <pcu_util.h>
#include "parma_centroids.h"
#include "parma_sides.h"

namespace {
  double getEntWeight(apf::Mesh* m, apf::MeshTag* w, apf::MeshEntity* e) {
    PCU_ALWAYS_ASSERT(m->hasTag(e,w));
    double entW = 0;
    m->getDoubleTag(e,w,&entW);
    return entW;
  }
  double selfWeight(apf::Mesh* m, apf::MeshTag* w) {
    const int dim = m->getDimension();
    apf::MeshIterator* it = m->begin(dim);
    apf::MeshEntity* e;
    double sum = 0;
    while ((e = m->iterate(it)))
      sum += getEntWeight(m,w,e);
    m->end(it);
    return sum;
  }
}

namespace parma {  
  Centroids::Centroids(apf::Mesh* m, apf::MeshTag* w, Sides* s) {
    weight = selfWeight(m,w);
    centroid = selfCentroid(m,w);
    init(m, s);
  }

  apf::Vector3 Centroids::self() {
    return centroid;
  }

  apf::Vector3 Centroids::selfCentroid(apf::Mesh* m, apf::MeshTag* w) {
    apf::Vector3 x(0,0,0);
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(m->getDimension());
    while ((e = m->iterate(it))) 
      x = x + (apf::getLinearCentroid(m, e) * getEntWeight(m,w,e));
    m->end(it);
    return x / weight;
  }

  void Centroids::init(apf::Mesh* m, Sides* s) {
    m->getPCU()->Begin();
    const Sides::Item* side;
    s->begin();
    while( (side = s->iterate()) ) 
      m->getPCU()->Pack(side->first, centroid);
    s->end();
    m->getPCU()->Send();
    while (m->getPCU()->Listen()) {
      apf::Vector3 otherCentroid;
      m->getPCU()->Unpack(otherCentroid);
      set(m->getPCU()->Sender(), otherCentroid);
    }
  }
} //end namespace
