#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <fstream>
#include <gmi.h>

#define MAX_ELEM_NODES 10

namespace apf {

static void ansys2apf(int ansysType, int& apfType, apf::FieldShape*& shape)
{
  switch (ansysType) {
    case 702:
      apfType = apf::Mesh::TET;
      shape = apf::getLagrange(1);
      return;
    case 902:
      apfType = apf::Mesh::TET;
      shape = apf::getLagrange(2);
      return;
  };
}

static void parseElemInt(std::string const& line, int at, int& out)
{
  assert(line.length() >= ((size_t)at + 1) * 6);
  std::string field = line.substr(at * 6, 6);
  if (field != "      ") {
    std::stringstream ss(field);
    ss >> out;
  }
}

static bool parseElem(std::istream& f, int nodes[MAX_ELEM_NODES],
    int& apfType, int& id, apf::FieldShape*& shape)
{
  std::string line;
  std::getline(f, line);
  if (!line.length())
    return false;
  for (int i = 0; i < 8; ++i)
    parseElemInt(line, i, nodes[i]);
  int ansysType;
  parseElemInt(line,  9, ansysType);
  parseElemInt(line, 13, id);
  ansys2apf(ansysType, apfType, shape);
  int nnodes = shape->getEntityShape(apfType)->countNodes();
  if (nnodes <= 8)
    return true;
  std::getline(f, line);
  for (int i = 0; i < (nnodes - 8); ++i)
    parseElemInt(line, i, nodes[i + 8]);
  return true;
}

typedef std::map<int, apf::Vector3> Nodes;

static bool parseNode(std::istream& f, int& id, apf::Vector3& p)
{
  std::string line;
  std::getline(f, line);
  if (!line.length())
    return false;
  std::stringstream ss(line);
  ss >> id;
  int i = 0;
  while (ss >> p[i])
    ++i;
  for (; i < 3; ++i)
    p[i] = 0;
  return true;
}

static void parseNodes(const char* nodefile, Nodes& nodes)
{
  std::ifstream f(nodefile);
  if (!f.is_open()) {
    fprintf(stderr, "couldn't open ANSYS node file \"%s\"\n", nodefile);
    abort();
  }
  std::pair<int, apf::Vector3> entry;
  while (parseNode(f, entry.first, entry.second))
    nodes.insert(entry);
}

typedef std::map<int, MeshEntity*> Vertices;

static Mesh2* parseElems(const char* elemfile, Nodes& nodes)
{
  Mesh2* m = 0;
  std::ifstream f(elemfile);
  if (!f.is_open()) {
    fprintf(stderr, "couldn't open ANSYS elem file \"%s\"\n", elemfile);
    abort();
  }
  Vertices verts;
  int en[MAX_ELEM_NODES];
  int type;
  int id;
  apf::FieldShape* shape = 0;
  apf::FieldShape* prevShape = 0;
  while (parseElem(f, en, type, id, shape)) {
    if (!m) {
      m = makeEmptyMdsMesh(gmi_load(".null"), Mesh::typeDimension[type], false);
      if (shape != m->getShape())
        changeMeshShape(m, shape, false);
    }
    if (prevShape)
      assert(prevShape == shape);
    int nev = Mesh::adjacentCount[type][0];
    int nen = shape->getEntityShape(type)->countNodes();
    Downward ev;
    for (int i = 0; i < nen; ++i)
      if (!nodes.count(en[i])) {
        fprintf(stderr, "node %d in file \"%s\" not found in node file\n",
            en[i], elemfile);
        abort();
      }
    int i = 0;
    for (i = 0; i < nev; ++i) {
      if (!verts.count(en[i])) {
        MeshEntity* v = m->createVert(0);
        m->setPoint(v, 0, nodes[en[i]]);
        verts[en[i]] = v;
      }
      ev[i] = verts[en[i]];
    }
    MeshEntity* e = buildElement(m, 0, type, ev);
    for (int d = 1; d <= m->getDimension(); ++d) {
      if (!shape->hasNodesIn(d))
        continue;
      Downward ee;
      int nee = m->getDownward(e, d, ee);
      for (int j = 0; j < nee; ++j) {
        assert(shape->countNodesOn(m->getType(ee[j])) == 1);
        m->setPoint(ee[j], 0, nodes[en[i]]);
        ++i;
      }
    }
    assert(i == nen);
    prevShape = shape;
  }
  return m;
}

Mesh2* loadMdsFromANSYS(const char* nodefile, const char* elemfile)
{
  Nodes nodes;
  parseNodes(nodefile, nodes);
  Mesh2* m = parseElems(elemfile, nodes);
  m->acceptChanges();
  deriveMdsModel(m);
  return m;
}

}
