#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfOmega_h.h>
#include <cstdlib>

#include <iostream>

#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  if (argc != 4) {
    if (PCU_Comm_Self() == 0) {
      std::cout << "\n";
      std::cout << "usage: osh2smb in.osh in.dmg out.smb\n";
      std::cout << "   or: osh2smb               (usage)\n";
    }
    PCU_Comm_Free();
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  gmi_register_null();
  gmi_model* model = gmi_load(argv[2]);
  {
    auto lib = Omega_h::Library(&argc, &argv);
    Omega_h::Mesh om(&lib);
    Omega_h::binary::read(argv[1], lib.world(), &om);
    apf::Mesh2* am = apf::makeEmptyMdsMesh(model, om.dim(), false);
    apf::from_omega_h(am, &om);
    am->writeNative(argv[3]);
    am->destroyNative();
    apf::destroyMesh(am);
  }
  PCU_Comm_Free();
  MPI_Finalize();
}
