#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfOmega_h.h>
#include <cstdlib>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  if (argc != 4) {
    if (PCU_Comm_Self() == 0) {
      std::cout << "\n";
      std::cout << "usage: smb2osh in.dmg in.smb out.osh\n";
      std::cout << "   or: smb2osh               (usage)\n";
    }
    PCU_Comm_Free();
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  gmi_register_null();
  apf::Mesh2* am = apf::loadMdsMesh(argv[1], argv[2]);
  {
    auto lib = Omega_h::Library(&argc, &argv);
    Omega_h::Mesh om(&lib);
    apf::to_omega_h(am, &om);
    am->destroyNative();
    apf::destroyMesh(am);
    Omega_h::binary::write(argv[3], &om);
  }
  PCU_Comm_Free();
  MPI_Finalize();
}
