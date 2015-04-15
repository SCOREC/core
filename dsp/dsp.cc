#include "dsp.h"
#include <apf.h>

namespace dsp {

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
    Boundary& fixed, Boundary& moving)
{
  smoother->smooth(df, fixed, moving);
  while (!tryToDisplace(m, df))
    adapter->adapt(m);
}

apf::Field* applyRigidMotion(apf::Mesh* m, Boundary& fixed,
    apf::Matrix3x3 const& r, apf::Vector3 const& t)
{
  apf::Field* dsp = apf::createFieldOn(m, "dsp", apf::VECTOR);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  while ((v = m->iterate(it)))
    if (fixed.count(m->toModel(v))) {
      apf::Vector3 x;
      m->getPoint(v, 0, x);
      apf::Vector3 nx;
      nx = (r * x) + t;
      apf::setVector(dsp, v, 0, nx - x);
    }
  m->end(it);
  return dsp;
}

}
