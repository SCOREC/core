#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_null.h>
#include <PCU.h>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_null();
  gmi_model* model = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 2, false);
  apf::Vector3 points[4] = {
    apf::Vector3(0,0,0),
    apf::Vector3(1,0,0),
    apf::Vector3(0,1,0),
    apf::Vector3(0,0,1)
  };
  apf::Vector3 param(0,0,0);
  apf::ModelEntity* surface = m->findModelEntity(2, 0);
  apf::MeshEntity* verts[4];
  for (int i = 0; i < 4; ++i)
    verts[i] = m->createVertex(surface, points[i], param);
  for (int i = 0; i < 4; ++i) {
    apf::MeshEntity* tv[3];
    for (int j = 0; j < 3; ++j)
      tv[j] = verts[apf::tet_tri_verts[i][j]];
    apf::buildElement(m, surface, apf::Mesh::TRIANGLE, tv);
  }
  m->acceptChanges();
  m->verify();
  apf::changeMdsDimension(m, 3);
  apf::ModelEntity* interior = m->findModelEntity(3, 0);
  apf::buildElement(m, interior, apf::Mesh::TET, verts);
  m->acceptChanges();
  m->verify();
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}


