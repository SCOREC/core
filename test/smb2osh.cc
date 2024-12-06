#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <lionPrint.h>
#include <apfOmega_h.h>
#include <cstdlib>

#include <iostream>

#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  if (argc != 4) {
    if (PCUObj.Self() == 0) {
      std::cout << "\n";
      std::cout << "usage: smb2osh in.dmg in.smb out.osh\n";
      std::cout << "   or: smb2osh               (usage)\n";
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  gmi_register_null();
  apf::Mesh2* am = apf::loadMdsMesh(argv[1], argv[2], &PCUObj);
  {
    auto lib = Omega_h::Library(&argc, &argv);
    Omega_h::Mesh om(&lib);
    apf::to_omega_h(&om, am);
    am->destroyNative();
    apf::destroyMesh(am);
    Omega_h::binary::write(argv[3], &om);
  }
  }
  MPI_Finalize();
}
