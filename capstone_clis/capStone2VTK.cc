#include <iostream>
#include <vector>

#include <PCU.h>
#include <pcu_util.h>
#include <gmi.h>
#include <gmi_cap.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfCAP.h>

#include "CapstoneModule.h"

using namespace CreateMG;
using namespace CreateMG::Attribution;
using namespace CreateMG::Mesh;
using namespace CreateMG::Geometry;


int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();

  if (argc != 3) {
    if(0==PCU_Comm_Self())
      std::cerr << "usage: " << argv[0]
        << " <cre file .cre> <output folder name>\n";
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  const char* creFileName = argv[1];
  const char* folderName = argv[2];

  // load capstone mesh
  CapstoneModule  cs("the_module", "Geometry Database : SMLIB",
    "Mesh Database : Create", "Attribution Database : Create");

  GeometryDatabaseInterface     *g = cs.get_geometry();
  MeshDatabaseInterface         *m = cs.get_mesh();

  PCU_ALWAYS_ASSERT(g);
  PCU_ALWAYS_ASSERT(m);

  M_GModel gmodel = cs.load_files(v_string(1, creFileName));

  M_MModel mmodel;
  std::vector<M_MModel> mmodels;
  MG_API_CALL(m, get_associated_mesh_models(gmodel, mmodels));
  std::cout << "Num Mesh Models loaded: " << mmodels.size() << std::endl;
  v_string names(mmodels.size());
  for (size_t i = 0; i < mmodels.size(); ++i) {
    MG_API_CALL(m, get_model_name(mmodels[i], names[i]));
    std::cout << i + 1 << ". " << names[i] << std::endl;
    std::string info;
    MG_API_CALL(m, print_info(mmodels[i], info));
    std::cout << info << std::endl;
  }
  int choice = 1;
  if (mmodels.size() > 1) {
    do {
      std::cout << "Choose 1-" << names.size() << ": ";
      std::cin >> choice;
    } while (!(1 <= choice && choice <= int(names.size())));
  }
  mmodel = mmodels[choice - 1];
  MG_API_CALL(m, set_active_model(mmodel));

  /* SET THE ADJACENCIES */
  MG_API_CALL(m, set_adjacency_state(REGION2FACE|
                                     REGION2EDGE|
                                     REGION2VERTEX|
                                     FACE2EDGE|
                                     FACE2VERTEX));
  MG_API_CALL(m, set_reverse_states());
  MG_API_CALL(m, set_adjacency_scope(TOPO_EDGE, SCOPE_FULL));
  MG_API_CALL(m, set_adjacency_scope(TOPO_FACE, SCOPE_FULL));
  MG_API_CALL(m, compute_adjacency());

  gmi_cap_start();
  gmi_register_cap();

  // convert the mesh to apf/mds mesh
  apf::Mesh2* mesh = apf::createMesh(m,g);

  apf::writeVtkFiles(folderName, mesh);

  gmi_cap_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}
