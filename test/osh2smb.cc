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
#ifndef SCOREC_NO_MPI
  MPI_Init(&argc, &argv);
#else
  (void) argc, (void) argv;
#endif
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  if (argc != 4) {
    if (PCUObj.Self() == 0) {
      std::cout << "\n";
      std::cout << "usage: osh2smb in.osh in.dmg out.smb\n";
      std::cout << "   or: osh2smb               (usage)\n";
    }
#ifndef SCOREC_NO_MPI
    MPI_Finalize();
#endif
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  gmi_register_null();
  gmi_model* model = gmi_load(argv[2]);
  {
    auto lib = Omega_h::Library(&argc, &argv);
    Omega_h::Mesh om(&lib);
    Omega_h::binary::read(argv[1], lib.world(), &om);
    apf::Mesh2* am = apf::makeEmptyMdsMesh(model, om.dim(), false, &PCUObj);
    apf::from_omega_h(am, &om);
    am->writeNative(argv[3]);
    am->destroyNative();
    apf::destroyMesh(am);
  }
  }
#ifndef SCOREC_NO_MPI
  MPI_Finalize();
#endif
}
