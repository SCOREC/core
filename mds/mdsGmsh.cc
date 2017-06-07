#include "apf.h"
#include "apfMDS.h"
#include "apfMesh2.h"
#include "apfShape.h"
#include "gmi.h" /* this is for gmi_getline... */

#include <cstdio>
#include <cstring>
#include <pcu_util.h>
#include <cstdlib>

namespace {

bool isQuadratic(int gmshType)
{
  switch (gmshType) {
    case 8: return true;
    case 9: return true;
    case 11: return true;
    default: return false;
  }
}

int apfFromGmsh(int gmshType)
{
  switch (gmshType) {
    case 1: return apf::Mesh::EDGE;
    case 2: return apf::Mesh::TRIANGLE;
    case 3: return apf::Mesh::QUAD;
    case 4: return apf::Mesh::TET;
    case 5: return apf::Mesh::HEX;
    case 6: return apf::Mesh::PRISM;
    case 7: return apf::Mesh::PYRAMID;
    case 8: return apf::Mesh::EDGE;
    case 9: return apf::Mesh::TRIANGLE;
    case 11: return apf::Mesh::TET;
    case 15: return apf::Mesh::VERTEX;
    default: return -1;
  }
}

struct Node {
  Node() { entity = 0; }
  apf::MeshEntity* entity;
  apf::Vector3 point;
};

struct Reader {
  apf::Mesh2* mesh;
  FILE* file;
  char* line;
  char* word;
  size_t linecap;
  bool isQuadratic;
  std::map<long, Node> nodeMap;
  std::map<long, apf::MeshEntity*> entMap[4];
};

void initReader(Reader* r, apf::Mesh2* m, const char* filename)
{
  r->mesh = m;
  r->file = fopen(filename, "r");
  if (!r->file) {
    fprintf(stderr,"couldn't open Gmsh file \"%s\"\n",filename);
    abort();
  }
  r->line = static_cast<char*>(malloc(1));
  r->line[0] = '\0';
  r->linecap = 1;
  r->isQuadratic = false;
}

void freeReader(Reader* r)
{
  free(r->line);
  fclose(r->file);
}

void getLine(Reader* r)
{
  ssize_t ret = gmi_getline(&r->line, &r->linecap, r->file);
  PCU_ALWAYS_ASSERT(ret != -1);
  r->word = r->line;
}

long getLong(Reader* r)
{
  long x;
  int nchars;
  int ret = sscanf(r->word, "%ld%n", &x, &nchars);
  PCU_ALWAYS_ASSERT(ret == 1);
  r->word += nchars;
  return x;
}

bool startsWith(char const* prefix, char const* s)
{
  int ls = strlen(s);
  int lp = strlen(prefix);
  if (ls < lp)
    return false;
  return strncmp(s, prefix, lp) == 0;
}

void seekMarker(Reader* r, char const* marker)
{
  while (!startsWith(marker, r->line))
    getLine(r);
  getLine(r);
}

void checkMarker(Reader* r, char const* marker)
{
  PCU_ALWAYS_ASSERT(startsWith(marker, r->line));
}

void readNode(Reader* r)
{
  long id;
  Node n;
  apf::Vector3& p = n.point;
  sscanf(r->line, "%ld %lf %lf %lf", &id, &p[0], &p[1], &p[2]);
  r->nodeMap[id] = n;
  getLine(r);
}

void readNodes(Reader* r)
{
  seekMarker(r, "$Nodes");
  long n = getLong(r);
  getLine(r);
  for (long i = 0; i < n; ++i)
    readNode(r);
  checkMarker(r, "$EndNodes");
}

apf::MeshEntity* lookupVert(Reader* r, long nodeId, apf::ModelEntity* g)
{
  PCU_ALWAYS_ASSERT(r->nodeMap.count(nodeId));
  Node& n = r->nodeMap[nodeId];
  if (n.entity)
    return n.entity;
  n.entity = r->mesh->createVert(g);
  r->mesh->setPoint(n.entity, 0, n.point);
  return n.entity;
}

void readElement(Reader* r)
{
  long id = getLong(r);
  long gmshType = getLong(r);
  if (isQuadratic(gmshType))
    r->isQuadratic = true;
  int apfType = apfFromGmsh(gmshType);
  PCU_ALWAYS_ASSERT(0 <= apfType);
  int nverts = apf::Mesh::adjacentCount[apfType][0];
  int dim = apf::Mesh::typeDimension[apfType];
  long ntags = getLong(r);
  PCU_ALWAYS_ASSERT(ntags >= 2);
  getLong(r); /* discard physical type */
  long gtag = getLong(r);
  for (long i = 2; i < ntags; ++i)
    getLong(r); /* discard all other element tags */
  apf::ModelEntity* g = r->mesh->findModelEntity(dim, gtag);
  apf::Downward verts;
  for (int i = 0; i < nverts; ++i) {
    long nid = getLong(r);
    verts[i] = lookupVert(r, nid, g);
  }
  if (dim != 0) {
    if (dim > r->mesh->getDimension())
      apf::changeMdsDimension(r->mesh, dim);
    apf::MeshEntity* ent = apf::buildElement(r->mesh, g, apfType, verts);
    if (r->isQuadratic)
      r->entMap[dim][id] = ent;
  }
  getLine(r);
}

void readElements(Reader* r)
{
  seekMarker(r, "$Elements");
  long n = getLong(r);
  getLine(r);
  for (long i = 0; i < n; ++i)
    readElement(r);
  checkMarker(r, "$EndElements");
}

static const double gmshTet10EdgeIndices[6] = {0, 1, 2, 3, 5, 4};

static int getQuadGmshIdx(const int apfIdx, const int apfType) {
  switch (apfType) {
    case apf::Mesh::TET :
      return gmshTet10EdgeIndices[apfIdx];
      break;
    default:
      return apfIdx;
  }
}

void readQuadraticElement(Reader* r)
{
  long id = getLong(r);
  long gmshType = getLong(r);
  int apfType = apfFromGmsh(gmshType);
  if (apfType == apf::Mesh::VERTEX) return;
  PCU_ALWAYS_ASSERT(0 <= apfType);
  PCU_ALWAYS_ASSERT_VERBOSE(isQuadratic(gmshType),
      "no support for variable p-order meshes");
  int nverts = apf::Mesh::adjacentCount[apfType][0];
  int nedges = apf::Mesh::adjacentCount[apfType][1];
  int dim = apf::Mesh::typeDimension[apfType];
  long ntags = getLong(r);
  getLong(r); /* discard physical type */
  getLong(r); /* discard geometric tag */
  for (long i = 2; i < ntags; ++i)
    getLong(r); /* discard all other tags */
  for (long i = 0; i < nverts; ++i)
    getLong(r); /* discard all vertex nodes */
  apf::Downward edges;
  apf::MeshEntity* ent = r->entMap[dim][id];
  r->mesh->getDownward(ent, 1, edges);
  apf::Field* coord = r->mesh->getCoordinateField();
  std::vector<long> nids(nedges);
  for (long i = 0; i < nedges; ++i)
    nids[i] = getLong(r);
  for (long i = 0; i < nedges; ++i) {
    long nid = nids[ getQuadGmshIdx(i, apfType) ];
    apf::Vector3 point = r->nodeMap[nid].point;
    apf::setVector(coord, edges[i], 0, point);
  }
}

void readQuadratic(Reader* r, apf::Mesh2* m, const char* filename)
{
  m->changeShape(apf::getSerendipity());
  initReader(r, m, filename);
  seekMarker(r, "$Elements");
  long n = getLong(r);
  getLine(r);
  for (long i = 0; i < n; ++i) {
    readQuadraticElement(r);
    getLine(r);
  }
  checkMarker(r, "$EndElements");
  freeReader(r);
}

void readGmsh(apf::Mesh2* m, const char* filename)
{
  Reader r;
  initReader(&r, m, filename);
  readNodes(&r);
  readElements(&r);
  freeReader(&r);
  m->acceptChanges();
  if (r.isQuadratic)
    readQuadratic(&r, m, filename);
}

}

namespace apf {

Mesh2* loadMdsFromGmsh(gmi_model* g, const char* filename)
{
  Mesh2* m = makeEmptyMdsMesh(g, 0, false);
  readGmsh(m, filename);
  return m;
}

}
