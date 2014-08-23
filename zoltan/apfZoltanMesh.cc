/*
 * Copyright (C) 2014 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfZoltanMesh.h"
#include "apfZoltanCallbacks.h"
#include "apfZoltan.h"
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
  opposite = 0;
}

ZoltanMesh::~ZoltanMesh()
{
  if (local)
    destroyNumbering(local);
  if (global)
    destroyGlobalNumbering(global);
  if (opposite) {
    const int sideDim = mesh->getDimension() - 1;
    removeTagFromDimension(mesh, opposite, sideDim);
    mesh->destroyTag(opposite);
  }
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
  assert(i == b->elements.getSize());
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
  if (!isLocal)
    opposite = tagOpposites(global, "zb_opposite");
  ZoltanData ztn(this);
  ztn.run();
  return convertResult(this, &ztn);
}

}
