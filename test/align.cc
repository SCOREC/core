#include <apfMesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apf.h>

void testTriEdge()
{
  gmi_model* model = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 2, false);
  apf::MeshEntity* v[3];
  for (int i = 0; i < 3; ++i)
    v[i] = m->createVert(0);
  apf::MeshEntity* ev[2];
  ev[0] = v[2]; ev[1] = v[1];
  apf::MeshEntity* e =
    apf::buildElement(m, 0, apf::Mesh::EDGE, ev);
  apf::MeshEntity* t =
    apf::buildElement(m, 0, apf::Mesh::TRIANGLE, v);
  int which, rotate;
  bool flip;
  apf::getAlignment(m, t, e, which, flip, rotate);
  assert(which == 1);
  assert(flip == true);
  assert(rotate == 0);
  m->destroyNative();
  apf::destroyMesh(m);
}

void testTetTri()
{
  gmi_model* model = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 3, false);
  apf::MeshEntity* v[4];
  for (int i = 0; i < 4; ++i)
    v[i] = m->createVert(0);
  apf::MeshEntity* tv[3];
  tv[0] = v[2]; tv[1] = v[1]; tv[2] = v[0];
  apf::MeshEntity* e =
    apf::buildElement(m, 0, apf::Mesh::TRIANGLE, tv);
  apf::MeshEntity* t =
    apf::buildElement(m, 0, apf::Mesh::TET, v);
  int which, rotate;
  bool flip;
  apf::getAlignment(m, t, e, which, flip, rotate);
  assert(which == 0);
  assert(flip == true);
  assert(rotate == 0);
  m->destroyNative();
  apf::destroyMesh(m);
}

int main()
{
  MPI_Init(0,0);
  PCU_Comm_Init();
  gmi_register_null();
  testTriEdge();
  testTetTri();
  PCU_Comm_Free();
  MPI_Finalize();
}
