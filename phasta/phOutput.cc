#include "phOutput.h"
#include "phLinks.h"
#include "phAdjacent.h"
#include <PCU.h>

namespace ph {

static void getCounts(Output& o)
{
  for (int i = 0; i < 4; ++i)
    o.nGlobalEntities[i] = apf::countOwned(o.mesh, i);
  PCU_Add_Ints(o.nGlobalEntities, 4);
  o.nOwnedNodes = apf::countOwned(o.mesh, 0);
  o.nOverlapNodes = o.mesh->count(0);
  o.nGlobalNodes = o.nGlobalEntities[0];
}

static void getGlobal(Output& o)
{
  apf::Mesh* m = o.mesh;
  apf::Numbering* n = apf::numberOwnedNodes(o.mesh, "ph_owned");
  apf::GlobalNumbering* gn = apf::makeGlobal(n);
  apf::synchronize(gn);
  o.arrays.globalNodeNumbers = new int[m->count(0)];
  apf::MeshEntity* v;
  size_t i = 0;
  apf::MeshIterator* it = m->begin(0);
  while ((v = m->iterate(it)))
    o.arrays.globalNodeNumbers[i++] = apf::getNumber(gn, apf::Node(v, 0));
  m->end(it);
  assert(i == m->count(0));
  apf::destroyGlobalNumbering(gn);
}

static void getBlocks(Output& o, ModelBounds& modelFaces)
{
  getAllBlocks(o.mesh, o.blocks, modelFaces);
}

static void getLinks(Output& o, apf::Numbering* n)
{
  Links links;
  getVertexLinks(o.mesh, links);
  size_t size;
  encodeLinks(n, links, size, o.arrays.ilwork);
  o.nlwork = size;
}

static void getInterior(Output& o, apf::Numbering* n)
{
  apf::Mesh* m = o.mesh;
  Blocks& bs = o.blocks.interior;
  int*** ien = new int**[bs.getSize()];
  for (int i = 0; i < bs.getSize(); ++i) {
    ien[i] = new int*[bs.nElements[i]];
    int t = bs.keys[i].elementType;
    int nv = bs.keys[i].nElementVertices;
    apf::MeshEntity* e;
    int j = 0;
    apf::MeshIterator* it = m->begin(m->getDimension());
    while ((e = m->iterate(it))) {
      if (getPhastaType(m, e) != t)
        continue;
      ien[i][j] = new int[nv];
      apf::Downward v;
      getVertices(m, e, v);
      for (int k = 0; k < nv; ++k)
        ien[i][j][k] = getNumber(n, v[k], 0, 0);
      ++j;
    }
    m->end(it);
  }
  o.arrays.ien = ien;
}

static void getBoundary(Output& o, ModelBounds& modelFaces, apf::Numbering* n)
{
  apf::Mesh* m = o.mesh;
  Blocks& bs = o.blocks.boundary;
  int*** ienb = new int**[bs.getSize()];
  int i = 0;
  APF_ITERATE(ModelBounds, modelFaces, mit) {
    apf::ModelEntity* mf = *mit;
    ienb[i] = new int*[bs.nElements[i]];
    int t = bs.keys[i].elementType;
    int nv = bs.keys[i].nElementVertices;
    apf::MeshEntity* f;
    int j = 0;
    apf::MeshIterator* it = m->begin(m->getDimension() - 1);
    while ((f = m->iterate(it))) {
      if (m->toModel(f) != mf)
        continue;
      apf::MeshEntity* e = m->getUpward(f, 0);
      if (getPhastaType(m, e) != t)
        continue;
      ienb[i][j] = new int[nv];
      apf::Downward v;
      getBoundaryVertices(m, e, f, v);
      for (int k = 0; k < nv; ++k)
        ienb[i][j][k] = getNumber(n, v[k], 0, 0);
      ++j;
    }
    m->end(it);
    ++i;
  }
  o.arrays.ienb = ienb;
}

void generateOutput(Input& in, apf::Mesh* mesh, Output& o)
{
  ModelBounds modelFaces; //FIXME: BC application not done yet
  o.in = &in;
  o.mesh = mesh;
  getCounts(o);
  getGlobal(o);
  getBlocks(o, modelFaces);
  apf::Numbering* n = apf::numberOverlapNodes(mesh, "ph_local");
  getLinks(o, n);
  getInterior(o, n);
  getBoundary(o, modelFaces, n);
  apf::destroyNumbering(n);
}

}
