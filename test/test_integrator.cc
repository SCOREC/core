#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <pcu_util.h>
#include <mpi.h>
#include <iostream>
#include <gmi_mesh.h>
#include <gmi_null.h>

class CountIntegrator : public apf::Integrator {
  protected:
    unsigned int numEnt;
  public:
    unsigned int getCount() {return numEnt;}
    void resetCount() { numEnt = 0; }
    CountIntegrator() : Integrator(1), numEnt(0) {};
    void inElement(apf::MeshElement *) 
    {
      numEnt++;
    }
    void atPoint(apf::Vector3 const& , double , double ) {
    }
};
int main(int argc, char ** argv) {
  MPI_Init(&argc, &argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  // argument should be model, mesh
  PCU_ALWAYS_ASSERT(argc == 3);

  gmi_register_mesh();
  gmi_register_null();
  apf::Mesh2* mesh = apf::loadMdsMesh(argv[1], argv[2], &PCUObj);
  CountIntegrator * countInt = new CountIntegrator();
  // test integration over implicitly defined mesh dimension
  countInt->process(mesh);
  PCU_ALWAYS_ASSERT(mesh->count(3) == countInt->getCount());
  // cannot test over vertices because integrator segfaults
  // when using a vertex "element"
  // test integration over explicitly defined mesh dimension
  for(int i=3; i>0; --i) {
    countInt->resetCount();
    countInt->process(mesh, i);
    PCU_ALWAYS_ASSERT(mesh->count(i) == countInt->getCount());
  }

  delete countInt;
  mesh->destroyNative();
  apf::destroyMesh(mesh);
  }
  MPI_Finalize();
  return 0;
}
