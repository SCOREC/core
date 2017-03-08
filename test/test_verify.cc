#include <apf.h>
#include <apfShape.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <PCU.h>
#include <vector>
#include <sstream>
#include <pcu_util.h>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  PCU_ALWAYS_ASSERT(argc == 3);
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  apf::DynamicArray<apf::Field*> fields(10);
  apf::FieldShape* shapes[10] = {
    apf::getLagrange(1),
    apf::getSerendipity(),
    apf::getHierarchic(2),
    apf::getIPShape(3,1),
    apf::getIPShape(3,2),
    apf::getVoronoiShape(3,1),
    apf::getVoronoiShape(3,2),
    apf::getIPFitShape(3,1),
    apf::getIPFitShape(3,2),
    apf::getConstant(3) };
  for (size_t i=0; i < fields.size(); ++i) {
    std::ostringstream oss;
    oss << "field_" << i;
    std::string name = oss.str();
    fields[i] = apf::createField(m, name.c_str(), apf::SCALAR, shapes[i]);
    apf::zeroField(fields[i]);
  }
  apf::verify(m);
  m->changeShape(shapes[1]);
  apf::verify(m);
  for (size_t i=0; i < fields.size(); ++i)
    apf::destroyField(fields[i]);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
