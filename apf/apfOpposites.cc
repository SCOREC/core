/*
 * Copyright (C) 2025 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <apfMesh.h>
#include <apfNumbering.h>
#include <pcu_util.h>

namespace apf {

static bool hasOtherSide(Mesh* m, MeshEntity* s)
{
  if (m->isShared(s))
    return true;
  if (!m->hasMatching())
    return false;
  Matches matches;
  m->getMatches(s, matches);
  return matches.getSize();
}

static Copy getOtherSide(Mesh* m, MeshEntity* s)
{
  if (m->isShared(s))
    return getOtherCopy(m, s);
  Matches matches;
  m->getMatches(s, matches);
  PCU_ALWAYS_ASSERT(matches.getSize() == 1);
  return matches[0];
}

static MeshEntity* getSideElement(Mesh* m, MeshEntity* s)
{
  return m->getUpward(s, 0);
}

static long getElementGid(GlobalNumbering* gn, MeshEntity* e)
{
  return getNumber(gn, Node(e, 0));
}

static void packOtherGid(GlobalNumbering* gn, MeshEntity* s)
{
  Mesh* m = getMesh(gn);
  Copy other = getOtherSide(m, s);
  long gid = getElementGid(gn, getSideElement(m, s));
  m->getPCU()->Pack(other.peer, other.entity);
  m->getPCU()->Pack(other.peer, gid);
}

static void unpackOtherGid(Mesh* m, MeshTag* t)
{
  MeshEntity* s;
  m->getPCU()->Unpack(s);
  long gid;
  m->getPCU()->Unpack(gid);
  m->setLongTag(s, t, &gid);
}

MeshTag* tagOpposites(GlobalNumbering* gn, const char* name)
{
  Mesh* m = getMesh(gn);
  MeshTag* t = m->createLongTag(name, 1);
  int sideDim = m->getDimension() - 1;
  m->getPCU()->Begin();
  MeshIterator* it = m->begin(sideDim);
  MeshEntity* e;
  while ((e = m->iterate(it)))
    if (hasOtherSide(m, e))
      packOtherGid(gn, e);
  m->end(it);
  m->getPCU()->Send();
  while (m->getPCU()->Receive())
    unpackOtherGid(m, t);
  return t;
}

} // namespace apf
