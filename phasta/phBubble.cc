#include "phBubble.h"
#include "phInput.h"
#include <apfMesh.h>
#include <apf.h>

namespace ph {

struct Bubble {
  int id;
  apf::Vector3 center;
  double radius;
};

typedef std::vector<Bubble> Bubbles;

void readBubbles(Bubbles& bubbles)
{
  (void)bubbles; //silence clang warning while development is done
  /* open file, read content, etc... */
}

void setBubbleScalars(apf::Mesh* m, apf::MeshEntity* v,
    Bubbles& bubbles, double* sol)
{
  apf::Vector3 v_center;
  m->getPoint(v, 0, v_center);
  /* search through bubbles, etc... */
  sol[5] = 42;
  sol[6] = 42;
}

void initBubbles(apf::Mesh* m, Input& in)
{
  Bubbles bubbles;
  readBubbles(bubbles);
  assert(in.ensa_dof >= 7);
  apf::NewArray<double> s(in.ensa_dof);
  apf::Field* f = m->findField("solution");
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  while ((v = m->iterate(it))) {
    apf::getComponents(f, v, 0, &s[0]);
    setBubbleScalars(m, v, bubbles, &s[0]);
    apf::setComponents(f, v, 0, &s[0]);
  }
  m->end(it);
}

}

