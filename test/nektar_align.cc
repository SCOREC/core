#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apf.h>
#include <PCU.h>
#include <vector>
#include <algorithm>

namespace apf {
/* the more dangerous a function is,
 * the more it gets used.
 */
void hackMdsAdjacency(Mesh2* in, MeshEntity* up, int i, MeshEntity* down);
}

struct VertId {
  apf::MeshEntity* vert;
  long id;
  VertId(apf::MeshEntity* a, long b):
    vert(a),
    id(b)
  {
  }
  bool operator<(const VertId& other) const
  {
    return id < other.id;
  }
};

typedef std::vector<VertId> VertIdArray;

static void getVertIds(apf::GlobalNumbering* n, apf::MeshEntity* e, VertIdArray& vids)
{
  apf::Mesh* m = apf::getMesh(n);
  apf::Downward v;
  int nv = m->getDownward(e, 0, v);
  vids.clear();
  vids.reserve(nv);
  for (int i = 0; i < nv; ++i)
    vids.push_back(VertId(v[i], apf::getNumber(n,apf::Node(v[i],0))));
}

static void alignTet(apf::Mesh2* m, apf::GlobalNumbering* n, apf::MeshEntity* tet)
{
  apf::MeshEntity* f[4];
  m->getDownward(tet, 2, f);
  VertIdArray new_v;
  getVertIds(n, tet, new_v);
  std::sort(new_v.begin(), new_v.end());
  apf::Vector3 p[4];
  for (int i = 0; i < 4; ++i)
    m->getPoint(new_v[i].vert, 0, p[i]);
  if ((p[3] - p[0]) * apf::cross((p[1] - p[0]), (p[2] - p[0])) < 0)
    std::swap(new_v[0], new_v[1]);
  apf::MeshEntity* new_f[4];
  for (int i = 0; i < 4; ++i) {
    apf::MeshEntity* fv[3];
    for (int j = 0; j < 3; ++j)
      fv[j] = new_v[apf::tet_tri_verts[i][j]].vert;
    new_f[i] = apf::findElement(m, apf::Mesh::TRIANGLE, fv);
  }
  for (int i = 0; i < 4; ++i)
    apf::hackMdsAdjacency(m, tet, i, new_f[i]);
}

static void alignTri(apf::Mesh2* m, apf::GlobalNumbering* n, apf::MeshEntity* tri)
{
  apf::MeshEntity* e[3];
  m->getDownward(tri, 1, e);
  VertIdArray new_v;
  getVertIds(n, tri, new_v);
  std::sort(new_v.begin(), new_v.end());
  apf::MeshEntity* new_e[3];
  for (int i = 0; i < 3; ++i) {
    apf::MeshEntity* ev[2];
    for (int j = 0; j < 2; ++j)
      ev[j] = new_v[apf::tri_edge_verts[i][j]].vert;
    new_e[i] = apf::findElement(m, apf::Mesh::EDGE, ev);
  }
  for (int i = 0; i < 3; ++i)
    apf::hackMdsAdjacency(m, tri, i, new_e[i]);
}

static void alignForNektar(apf::Mesh2* m)
{
  apf::GlobalNumbering* n = apf::makeGlobal(apf::numberOwnedNodes(m, "nektar_id"));
  apf::synchronize(n);
  apf::MeshIterator* it;
  apf::MeshEntity* e;
  it = m->begin(2);
  while ((e = m->iterate(it)))
    alignTri(m, n, e);
  m->end(it);
  it = m->begin(3);
  while ((e = m->iterate(it)))
    alignTet(m, n, e);
  m->end(it);
}

int main(int argc, char** argv)
{
  assert(argc==4);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  alignForNektar(m);
  m->writeNative(argv[3]);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
