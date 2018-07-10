#include <PCU.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfAggregateNumbering.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <pcu_util.h>

void test_numbering(apf::Mesh * m)
{
  apf::Field * f1 = apf::createPackedField(m, "u1", 3);
  apf::Field * f2 = apf::createPackedField(m, "u2", 3);
  apf::Field * f3 = apf::createPackedField(m, "u3", 3);
  apf::Field * f4 = apf::createPackedField(m, "u4", 3);
  apf::Sharing * def = apf::getSharing(m);
  // global scope natural order
  apf::AggNumbering * num_gbl_nat = apf::createAggNumbering(f1,1,3,MPI_COMM_WORLD,def);
  // global scope block order
  apf::AggNumbering * num_gbl_blk = apf::createAggNumbering(f2,3,1,MPI_COMM_WORLD,def);
  // local scope natural order
  apf::AggNumbering * num_lcl_nat = apf::createAggNumbering(f3,1,3,MPI_COMM_SELF);
  // local scope block order
  apf::AggNumbering * num_lcl_blk = apf::createAggNumbering(f4,3,1,MPI_COMM_SELF);
}

void write_output(apf::Mesh * m, const char * fn)
{
  apf::writeVtkFiles(fn,m);
}

int main(int ac, char * av[])
{
  PCU_ALWAYS_ASSERT(ac == 4);
  MPI_Init(&ac,&av);
  PCU_Comm_Init();
  gmi_sim_start();
  gmi_register_sim();
  gmi_register_mesh();
  apf::Mesh2 * m = apf::loadMdsMesh(av[1],av[2]);
  apf::reorderMdsMesh(m);
  test_numbering(m);
  write_output(m,av[3]);
  m->destroyNative();
  apf::destroyMesh(m);
  gmi_sim_stop();
  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
}
