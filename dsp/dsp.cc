#include "dsp.h"
#include <apf.h>
#include <gmi.h>
#include <vector>

using namespace std;

namespace dsp {
  
  static void closeBoundaryRec(gmi_model* gm, gmi_ent* e, Boundary& b)
  {
    b.insert((apf::ModelEntity*)e);
    int dim = gmi_dim(gm, e);
    if (!dim)
      return;
    gmi_set* s = gmi_adjacent(gm, e, gmi_dim(gm, e) - 1);
    for (int i = 0; i < s->n; ++i)
      closeBoundaryRec(gm, s->e[i], b);
    gmi_free_set(s);
  }
  
  void closeBoundary(apf::Mesh* m, Boundary& b)
  {
    gmi_model* gm = m->getModel();
    Boundary nb;
    APF_ITERATE(Boundary, b, it)
    closeBoundaryRec(gm, (gmi_ent*) *it, nb);
    b = nb;
  }
  
  bool tryToDisplace(apf::Mesh2* m, apf::Field* df)
  {
    double f = 1;
    apf::axpy(f, df, m->getCoordinateField());
    if (0 == verifyVolumes(m, false))
      return true;
    do {
      f /= 2;
      apf::axpy(-f, df, m->getCoordinateField());
    } while (0 != verifyVolumes(m, false));
    apf::axpy(-f, df, df);
    return false;
  }
  
  void displace(apf::Mesh2* m, apf::Field* df,
                Smoother* smoother, Adapter* adapter,
                Boundary& fixed, Boundary& moving,
                vector < apf::MeshEntity* >& V_total,
                int& in_0, int& fb_0)
  {
    smoother->preprocess(m, fixed, moving, V_total, in_0, fb_0);
    smoother->smooth(df, V_total, in_0, fb_0);
    smoother->cleanup();
    while (!tryToDisplace(m, df))
      adapter->adapt(m);
  }
  
  apf::Field* applyRigidMotion(apf::Mesh* m, Boundary& moving,
                               apf::Matrix3x3 const& r, apf::Vector3 const& t)
  {
    apf::Field* dsp = apf::createFieldOn(m, "dsp", apf::VECTOR);
    apf::MeshIterator* it = m->begin(0);
    apf::MeshEntity* v;
    while ((v = m->iterate(it))) {
      if (moving.count(m->toModel(v))) {
        apf::Vector3 x;
        m->getPoint(v, 0, x);
        apf::Vector3 nx;
        nx = (r * x) + t;
        apf::setVector(dsp, v, 0, nx - x);
      } else {
        apf::setVector(dsp, v, 0, apf::Vector3(0,0,0));
      }
    }
    m->end(it);
    return dsp;
  }
  
}