#include "phOutput.h"
#include "phLinks.h"
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

static void getLinks(Output& o, apf::Numbering* n)
{
  Links links;
  getVertexLinks(o.mesh, links);
  size_t size;
  encodeLinks(n, links, size, o.arrays.ilwork);
  o.nlwork = size;
}

static void getGlobal(Output& o)
{
  apf::Mesh* m = o.mesh;
  apf::Numbering* n = apf::numberOwnedNodes(o.mesh, "ph_owned");
  apf::GlobalNumbering* gn = apf::makeGlobal(n);
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

static void getBlocks(Output& o)
{
  ModelBounds modelFaces; //FIXME: BC application not done yet
  getAllBlocks(o.mesh, o.blocks, modelFaces);
}

void generateOutput(Input& in, apf::Mesh* mesh, Output& o)
{
  o.in = &in;
  o.mesh = mesh;
  getCounts(o);
  getGlobal(o);
  getBlocks(o);
  apf::Numbering* n = apf::numberOverlapNodes(mesh, "ph_local");
  getLinks(o, n);
  apf::destroyNumbering(n);
}

}
