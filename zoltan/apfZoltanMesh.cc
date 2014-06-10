/*
 * Copyright (C) 2014 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfZoltanMesh.h"
#include "apfZoltanCallbacks.h"
#include <PCU.h>

namespace apf {

ZoltanMesh::ZoltanMesh(Mesh* mesh_, bool isLocal_, int method_, int approach_,
    bool dbg)
{
  mesh = mesh_;
  isLocal = isLocal_;
  method = method_;
  approach = approach_;
  debug = dbg;
  local = 0;
  global = 0;
  opposite = mesh->createLongTag("zb_opposite",1);
}

ZoltanMesh::~ZoltanMesh()
{
  if (local)
    destroyNumbering(local);
  if (global)
    destroyGlobalNumbering(global);
  const int sideDim = mesh->getDimension() - 1;
  removeTagFromDimension(mesh, opposite, sideDim);
  mesh->destroyTag(opposite);
}

static void getElements(ZoltanMesh* b)
{
  Mesh* m = b->mesh;
  b->elements.setSize(m->count(m->getDimension()));
  MeshIterator* it = m->begin(m->getDimension());
  MeshEntity* e;
  size_t i = 0;
  while ((e = m->iterate(it)))
    b->elements[i++] = e;
  assert(i = b->elements.getSize());
  m->end(it);
}

static void setupNumberings(ZoltanMesh* b)
{
  b->local = numberElements(b->mesh, "zoltan_element");
  if (!b->isLocal) {
    Numbering* tmp = numberElements(b->mesh, "zoltan");
    b->global = makeGlobal(tmp);
  }
}

static bool hasOther(Mesh* m, MeshEntity* s)
{
  if (m->isShared(s))
    return true;
  if (!m->hasMatching())
    return false;
  Matches matches;
  m->getMatches(s, matches);
  return matches.getSize();
}

static std::pair<int, MeshEntity*> getOther(Mesh* m, MeshEntity* s)
{
  if (m->isShared(s))
    return getOtherCopy(m, s);
  Matches matches;
  m->getMatches(s, matches);
  assert(matches.getSize() == 1);
  return std::make_pair(matches[0].peer, matches[0].entity);
}

static MeshEntity* getSideElement(Mesh* m, MeshEntity* s)
{
  return m->getUpward(s, 0);
}

static long getElementGid(ZoltanMesh* b, MeshEntity* e)
{
  return getNumber(b->global, Node(e, 0));
}

static void packOpposite(ZoltanMesh* b, MeshEntity* s)
{
  std::pair<int, MeshEntity*> other = getOther(b->mesh, s);
  long gid = getElementGid(b, getSideElement(b->mesh, s));
  PCU_COMM_PACK(other.first, other.second);
  PCU_COMM_PACK(other.first, gid);
}

static void unpackOpposite(ZoltanMesh* b)
{
  MeshEntity* s;
  PCU_COMM_UNPACK(s);
  long gid;
  PCU_COMM_UNPACK(gid);
  b->mesh->setLongTag(s, b->opposite, &gid);
}

static void tagOpposites(ZoltanMesh* b)
{
  
  if (b->isLocal)
    return;
  Mesh* m = b->mesh;
  int sideDim = m->getDimension() - 1;
  PCU_Comm_Begin();
  MeshIterator* it = m->begin(sideDim);
  MeshEntity* e;
  while ((e = m->iterate(it)))
    if (hasOther(m, e))
      packOpposite(b, e);
  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while (!PCU_Comm_Unpacked())
      unpackOpposite(b);
  m->end(it);
}

static Migration* convertResult(ZoltanMesh* b, ZoltanData* ztn)
{
  Migration* plan = new Migration(b->mesh);
  //fill out the plan from zoltan class ztn->get(int localId)
  for (int ind=0;ind<ztn->getNumExported();ind++) {
    int lid;
    int exportPart;
    ztn->getExport(ind,&lid,&exportPart);
    plan->send(b->elements[lid],exportPart);
  }
  return plan;
}

Migration* ZoltanMesh::run(MeshTag* w, double tol, int mult)
{
  weights = w;
  tolerance = tol;
  multiple = mult;
  setupNumberings(this);
  getElements(this);
  tagOpposites(this);
  ZoltanData ztn(this);
  ztn.run();
  return convertResult(this, &ztn);
}

}
