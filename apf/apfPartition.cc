#include "apfPartition.h"
#include "apfMesh2.h"
#include "apf.h"

namespace apf {

static void remapResidence(apf::Mesh2* m, Remap& remap)
{
  for (int d = 0; d <= m->getDimension(); ++d) {
    apf::MeshIterator* it = m->begin(d);
    apf::MeshEntity* e;
    while ((e = m->iterate(it))) {
      apf::Parts residence;
      m->getResidence(e, residence);
      apf::Parts newResidence;
      APF_ITERATE(apf::Parts, residence, rit)
        newResidence.insert( remap(*rit) );
      m->setResidence(e, newResidence);
    }
    m->end(it);
  }
}

static void remapRemotes(apf::Mesh2* m, Remap& remap)
{
  for (int d = 0; d < m->getDimension(); ++d) {
    apf::MeshIterator* it = m->begin(d);
    apf::MeshEntity* e;
    while ((e = m->iterate(it))) {
      if ( ! m->isShared(e))
        continue;
      apf::Copies remotes;
      m->getRemotes(e, remotes);
      apf::Copies newRemotes;
      APF_ITERATE(apf::Copies, remotes, rit)
        newRemotes[ remap(rit->first) ] = rit->second;
      m->setRemotes(e, newRemotes);
    }
    m->end(it);
  }
}

static void remapMatches(apf::Mesh2* m, Remap& remap)
{
  if (!m->hasMatching())
    return;
  for (int d = 0; d < m->getDimension(); ++d) {
    apf::MeshIterator* it = m->begin(d);
    apf::MeshEntity* e;
    while ((e = m->iterate(it))) {
      apf::Matches matches;
      m->getMatches(e, matches);
      if (!matches.getSize())
        continue;
      m->clearMatches(e);
      for (size_t i = 0; i < matches.getSize(); ++i)
        m->addMatch(e, remap( matches[i].peer ), matches[i].entity);
    }
    m->end(it);
  }
}

void remapPartition(apf::Mesh2* m, Remap& remap)
{
  remapResidence(m, remap);
  remapRemotes(m, remap);
  remapMatches(m, remap);
  m->acceptChanges();
}

}
