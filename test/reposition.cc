#include <gmi_null.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apf.h>
#include <PCU.h>
#include <maReposition.h>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
#if 0
  gmi_register_null();
  gmi_model* model = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 3, false);
  apf::Vector3 vx[4] =
  {apf::Vector3(0,0,0),
   apf::Vector3(1,0,0),
   apf::Vector3(0,1,0),
   apf::Vector3(0,0,1)};
  apf::MeshEntity* tet =
      apf::buildOneElement(m, m->findModelEntity(3, 0), apf::Mesh::TET, vx);
  apf::MeshEntity* tv[4];
  m->getDownward(tet, 0, tv);
  apf::MeshEntity* v = tv[0];
#else
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1], argv[2]);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  while ((v = m->iterate(it))) {
    if (m->getModelType(m->toModel(v)) == 3)
      break;
  }
  m->end(it);
  apf::Vector3 randx;
  for (int i = 0; i < 3; ++i)
    randx[i] = (rand() / ((double) RAND_MAX)) * 70 - 1;
  m->setPoint(v, 0, randx);
#endif
  apf::writeVtkFiles("before", m);
  ma::repositionVertex(m, v, 1.0);
  apf::writeVtkFiles("after", m);
  PCU_Comm_Free();
  MPI_Finalize();
}
