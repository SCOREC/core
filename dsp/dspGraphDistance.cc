#include "dspGraphDistance.h"
#include <pcu_util.h>

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
    }
  }
  m->end(it);
  for (int layer = 0; m->getPCU()->Or(!empty(vs, first)); ++layer) {
    size_t layer_end = vs.size();
    while (first < layer_end) {
      v = pop(vs, first);
      int prev_val = apf::getNumber(dst, v, 0, 0);
      PCU_ALWAYS_ASSERT(prev_val == layer);
      apf::Up up;
      m->getUp(v, up);
      apf::MeshEntity* ov;
      for (int i = 0; i < up.n; ++i) {
        ov = apf::getEdgeVertOppositeVert(m, up.e[i], v);
        if (!apf::isNumbered(dst, ov, 0, 0)) {
          apf::number(dst, ov, 0, 0, layer + 1);
          push(vs, ov);
        }
      }
    }
    m->getPCU()->Begin();
    apf::MeshEntity* sv;
    for (size_t i = first; i < vs.size(); ++i) {
      sv = vs[i];
      if (!m->isShared(sv))
        continue;
      int val = apf::getNumber(dst, sv, 0, 0);
      PCU_ALWAYS_ASSERT(val == layer + 1);
      apf::Copies remotes;
      m->getRemotes(sv, remotes);
      APF_ITERATE(apf::Copies, remotes, rit)
        m->getPCU()->Pack(rit->first, rit->second);
    }
    m->getPCU()->Send();
    while (m->getPCU()->Receive()) {
      m->getPCU()->Unpack(sv);
      if (!apf::isNumbered(dst, sv, 0, 0)) {
        apf::number(dst, sv, 0, 0, layer + 1);
        push(vs, sv);
      }
    }
  }
  PCU_ALWAYS_ASSERT(vs.size() == m->count(0));
  return dst;
}

}
