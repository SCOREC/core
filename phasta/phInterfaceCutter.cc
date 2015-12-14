#include "phInterfaceCutter.h"
#include <apfMDS.h>
#include <apf.h>
#include <ph.h>
#include <cassert>
#include <iostream>

namespace apf {
/* so dangerous is this function that it
   was never listed in a public header.
   don't call it unless you understand everything
   about meshes */
void hackMdsAdjacency(Mesh2* in, MeshEntity* up, int i, MeshEntity* down);
}

namespace ph {

typedef std::map<gmi_ent*, int> MaterialMap;
typedef std::set<int> MaterialSet;

static bool isInterface(gmi_model* gm, gmi_ent* ge, FieldBCs& fbcs)
{
  int d = gmi_dim(gm, ge);
  if (d > 2)
    return false;
  if (d == 2)
    return getBCValue(gm, fbcs, ge) != 0;
  bool out = false;
  gmi_set* s = gmi_adjacent(gm, ge, d + 1);
  for (int i = 0; i < s->n; ++i)
    if (isInterface(gm, s->e[i], fbcs)) {
      out = true;
      break;
    }
  gmi_free_set(s);
  return out;
}

/* depth first search marking with mid,
   geometric regions are graph nodes and
   geometric faces that are *not* interface
   faces are graph edges */
static void findMaterialsDFS(gmi_model* gm, FieldBCs& fbcs,
    MaterialMap& mm, gmi_ent* ge, int mid)
{
  mm[ge] = mid;
  gmi_set* fs = gmi_adjacent(gm, ge, 2);
  for (int i = 0; i < fs->n; ++i) {
    if (isInterface(gm, fs->e[i], fbcs))
      continue;
    gmi_set* rs = gmi_adjacent(gm, fs->e[i], 3);
    for (int j = 0; j < rs->n; ++j) {
      if (mm.count(rs->e[j]))
        continue;
      findMaterialsDFS(gm, fbcs, mm, rs->e[j], mid);
    }
    gmi_free_set(rs);
  }
  gmi_free_set(fs);
}

/* find connected components of model regions
   connected by non-interface faces,
   output is in the material map */
static void findMaterials(gmi_model* gm, FieldBCs& fbcs,
    MaterialMap& mm)
{
  int mid = 0;
  gmi_iter* rit = gmi_begin(gm, 3);
  gmi_ent* ge;
  while ((ge = gmi_next(gm, rit))) {
    if (mm.count(ge))
      continue;
    findMaterialsDFS(gm, fbcs, mm, ge, mid);
    ++mid;
  }
  gmi_end(gm, rit);
}

static void replaceAdjacencies(apf::Mesh2* m,
    apf::MeshEntity* of, apf::MeshEntity* olda, apf::MeshEntity* newa)
{
  int ad = apf::getDimension(m, olda);
  int od = apf::getDimension(m, of);
  if (od <= ad)
    return;
  apf::Downward down;
  int ndown = m->getDownward(of, od - 1, down);
  for (int i = 0; i < ndown; ++i) {
    if (down[i] == olda)
      hackMdsAdjacency(m, of, i, newa);
    else
      replaceAdjacencies(m, down[i], olda, newa);
  }
}

static apf::MeshEntity* cloneEntity(apf::Mesh2* m, apf::MeshEntity* e)
{
  int et = m->getType(e);
  int ed = apf::Mesh::typeDimension[et];
  apf::ModelEntity* me = m->toModel(e);
  if (ed == 0) {
    apf::Vector3 pt;
    apf::Vector3 pm;
    m->getPoint(e, 0, pt);
    m->getParam(e, pm);
    return m->createVertex(me, pt, pm);
  } else {
    apf::Downward down;
    m->getDownward(e, ed - 1, down);
    return m->createEntity(et, me, down);
  }
}

static void cutEntity(apf::Mesh2* m, MaterialMap& mm, apf::MeshEntity* e)
{
  apf::Adjacent elements;
  m->getAdjacent(e, m->getDimension(), elements);
  MaterialSet ms;
  APF_ITERATE(apf::Adjacent, elements, it)
    ms.insert(mm[ (gmi_ent*) (m->toModel(*it)) ]);
  std::vector<apf::MeshEntity*> ents;
  ents.reserve(ms.size());
  ms.erase(ms.begin());
  assert(ms.size());
  ents.push_back(e);
  APF_ITERATE(MaterialSet, ms, it) {
    apf::MeshEntity* ne = cloneEntity(m, e);
    APF_ITERATE(apf::Adjacent, elements, eit)
      if (mm[ (gmi_ent*) (m->toModel(*eit))] == *it)
        replaceAdjacencies(m, *eit, e, ne);
    ents.push_back(ne);
  }
  APF_ITERATE(std::vector<apf::MeshEntity*>, ents, ait)
    APF_ITERATE(std::vector<apf::MeshEntity*>, ents, bit)
      if (*ait != *bit)
        m->addMatch(*ait, m->getId(), *bit);
}

static void cutEntities(apf::Mesh2* m, FieldBCs& fbcs, MaterialMap& mm)
{
  apf::setMdsMatching(m, true);
  for (int d = m->getDimension() - 1; d >= 0; --d) {
    apf::MeshIterator* it = m->begin(d);
    apf::MeshEntity* e;
    std::vector<apf::MeshEntity*> toCut;
    toCut.reserve(m->count(d));
    while ((e = m->iterate(it))) {
      if (isInterface(m->getModel(), (gmi_ent*) m->toModel(e), fbcs))
        toCut.push_back(e);
    }
    m->end(it);
    for (size_t i = 0; i < toCut.size(); ++i)
      cutEntity(m, mm, toCut[i]);
  }
}

void cutInterface(apf::Mesh2* m, BCs& bcs)
{
  std::string name("DG interface");
  if (!haveBC(bcs, name))
    ph::fail("no DG interface attributes!");
  FieldBCs& fbcs = bcs.fields[name];
  MaterialMap mm;
  findMaterials(m->getModel(), fbcs, mm);
  cutEntities(m, fbcs, mm);
}

}
