#include "dspGraphDistance.h"
#include <PCU.h>

namespace dsp {

static void push(std::vector<apf::MeshEntity*>& vs,
    apf::MeshEntity* v)
{
  vs.push_back(v);
}

static bool empty(std::vector<apf::MeshEntity*> const& vs, size_t first)
{
  return first == vs.size();
}

static apf::MeshEntity* pop(std::vector<apf::MeshEntity*> const& vs,
    size_t& first)
{
  apf::MeshEntity* v = vs[first];
  ++first;
  return v;
}

apf::Numbering* getGraphDistance(apf::Mesh* m, Boundary& seed,
    std::vector<apf::MeshEntity*>& vs)
{
  apf::Numbering* dst = apf::createNumbering(m, "graph_dist", m->getShape(), 1);
  apf::MeshIterator* it;
  apf::MeshEntity* v;
  size_t first = 0;
  vs.reserve(m->count(0));
  it = m->begin(0);
  while ((v = m->iterate(it))) {
    if (seed.count(m->toModel(v))) {
      apf::number(dst, v, 0, 0, 0);
      push(vs, v);
    } else {
      apf::number(dst, v, 0, 0, -1);
    }
  }
  m->end(it);
  int layer = 0;
  while (PCU_Or(!empty(vs, first))) {
    v = pop(vs, first);
    int prev_val = apf::getNumber(dst, v, 0, 0);
    if (prev_val != layer) {
      PCU_Comm_Begin();
      apf::MeshEntity* sv;
      for (size_t i = first; i < vs.size(); ++i) {
        sv = vs[i];
        if (!m->isShared(sv))
          continue;
        int val = apf::getNumber(dst, sv, 0, 0);
        assert(val == layer + 1);
        apf::Copies remotes;
        m->getRemotes(sv, remotes);
        APF_ITERATE(apf::Copies, remotes, rit)
          PCU_COMM_PACK(rit->first, rit->second);
      }
      PCU_Comm_Send();
      while (PCU_Comm_Receive()) {
        PCU_COMM_UNPACK(sv);
        int val = apf::getNumber(dst, sv, 0, 0);
        if (val == -1) {
          apf::number(dst, sv, 0, 0, layer + 1);
          push(vs, sv);
        }
      }
      ++layer;
    }
    assert(prev_val == layer);
    apf::Up up;
    m->getUp(v, up);
    apf::MeshEntity* ov;
    for (int i = 0; i < up.n; ++i) {
      ov = apf::getEdgeVertOppositeVert(m, up.e[i], v);
      int val = apf::getNumber(dst, ov, 0, 0);
      if (val == -1) {
        apf::number(dst, ov, 0, 0, layer + 1);
        push(vs, ov);
      }
    }
  }
  assert(vs.size() == m->count(0));
  return dst;
}

};
