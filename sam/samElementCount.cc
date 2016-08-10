#include "samElementCount.h"

#include <apf.h>
#include <apfMesh.h>
#include <PCU.h>

namespace sam {

/* The algorithms here are based on Section 2.7 of:
Pain, C. C., et al.
"Tetrahedral mesh optimisation and adaptivity for
 steady-state and transient finite element calculations."
Computer Methods in Applied Mechanics and Engineering
190.29 (2001): 3771-3796.

We make slight modifications because we do
not store the metric tensor directly, rather
we need to scale desired edge lengths, so we
use the inverse of the square root of the $\beta$ factor.
(a metric tensor eigenvalue is $(1/h^2)$, where $h$ is desired length).
*/

double getVolumeChange(int dim, double h) {
  if (dim == 2) return 1. / (h * h);
  if (dim == 3) return 1. / (h * h * h);
  apf::fail("bad dim");
}

class TotalMetricVolumeIso : public apf::Integrator {
  apf::Field* iso_field;
  int dim;
  apf::Element* element;

public:
  double sum;

  TotalMetricVolumeIso(apf::Field* iso_field_):
    apf::Integrator(1),
    iso_field(iso_field_),
    dim(apf::getMesh(iso_field_)->getDimension()),
    sum(0) {
  }
  virtual void inElement(apf::MeshElement* me) {
    element = apf::createElement(iso_field, me);
  }
  virtual void outElement() {
    apf::destroyElement(element);
  }
  virtual void atPoint(apf::Vector3 const& xi, double w, double dV) {
    double h = apf::getScalar(element, xi);
    double vhat = getVolumeChange(dim, h);
    sum += vhat * w * dV;
  }
  virtual void parallelReduce() {
    sum = PCU_Add_Double(sum);
  }
};

double getTotalMetricVolumeIso(apf::Field* iso_field) {
  apf::Mesh* m = apf::getMesh(iso_field);
  TotalMetricVolumeIso integrator(iso_field);
  integrator.process(m);
  return integrator.sum;
}

double getPerfectVolume(int dim) {
  if (dim == 2) return sqrt(3.) / 4.;
  if (dim == 3) return sqrt(2.) / 12.;
  apf::fail("bad dim");
}

double getVolumeScalar(int dim, double targetElementCount,
    double currentMetricVolume) {
  return getPerfectVolume(dim) * targetElementCount / currentMetricVolume;
}

double getLengthScalar(int dim, double volumeScalar) {
  if (dim == 2) return 1. / sqrt(volumeScalar);
  if (dim == 3) return 1. / pow(volumeScalar, 1. / 3.);
  apf::fail("bad dim");
}

double getLengthScalar(int dim, double targetElementCount,
    double currentMetricVolume) {
  return getLengthScalar(dim,
      getVolumeScalar(dim, targetElementCount, currentMetricVolume));
}

double getIsoLengthScalar(apf::Field* iso_field, double targetElementCount) {
  apf::Mesh* m = apf::getMesh(iso_field);
  double currentMetricVolume = getTotalMetricVolumeIso(iso_field);
  return getLengthScalar(m->getDimension(), targetElementCount,
      currentMetricVolume);
}

void scaleIsoSizeField(apf::Field* iso_field, double targetElementCount) {
  apf::Mesh* m = apf::getMesh(iso_field);
  double lengthScalar = getIsoLengthScalar(iso_field, targetElementCount);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* vert;
  while ((vert = m->iterate(it))) {
    apf::setScalar(iso_field, vert, 0,
        apf::getScalar(iso_field, vert, 0) * lengthScalar);
  }
  m->end(it);
}

}
