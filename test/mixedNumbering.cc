#include <PCU.h>
#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <apfNumbering.h>
#include <apfMixedNumbering.h>
#include <gmi_mesh.h>
#include <pcu_util.h>
#include <sstream>

static void test_numbering(apf::Mesh* m) {
  apf::FieldShape* S2 = apf::getSerendipity();
  apf::FieldShape* S1 = apf::getLagrange(1);
  apf::Field* f1 = apf::createField(m, "u", apf::VECTOR, S2);
  apf::Field* f2 = apf::createField(m, "p", apf::SCALAR, S1);
  apf::zeroField(f1);
  apf::zeroField(f2);
  std::vector<apf::Field*> fields;
  fields.push_back(f1);
  fields.push_back(f2);
  std::vector<apf::Numbering*> owned;
  std::vector<apf::Numbering*> ghost;
  std::vector<apf::GlobalNumbering*> global;
  int num_owned = apf::numberOwned(fields, owned);
  int num_ghost = apf::numberGhost(fields, ghost);
  apf::makeGlobal(owned, global);
  for (size_t n=0; n < owned.size(); ++n)
    apf::destroyNumbering(owned[n]);
  for (size_t n=0; n < global.size(); ++n)
    apf::synchronize(global[n]);
  PCU_Debug_Open();
  PCU_Debug_Print("number owned: %d\n", num_owned);
  PCU_Debug_Print("number ghost: %d\n", num_ghost);
}

static void write_output(apf::Mesh* m, const char* out) {
  std::ostringstream name1;
  name1 << out << "_linear";
  apf::writeVtkFiles(name1.str().c_str(), m);
  m->changeShape(apf::getSerendipity());
  std::ostringstream name2;
  name2 << out << "_quadratic";
  apf::writeVtkFiles(name2.str().c_str(), m);
}

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc==4);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  apf::reorderMdsMesh(m);
  test_numbering(m);
  write_output(m, argv[3]);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
