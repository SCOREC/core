#include <cstring>
#include <cstdlib>

// Output
#include <lionPrint.h>

// Parallelism
#include <PCU.h>
#include <pcu_util.h>

// Mesh interfaces
#include <apf.h>
#include <apfCAP.h>
#include <apfMDS.h>
#include <apfConvert.h>

// Geometry interfaces
#include <gmi.h>
#include <gmi_cap.h>

// Mesh adapt
#include <ma.h>

using namespace CreateMG;
using namespace CreateMG::Attribution;
using namespace CreateMG::Mesh;
using namespace CreateMG::Geometry;

#include "capVolSizeFields.h"

namespace {

void myExit(int exit_code = EXIT_SUCCESS) {
  gmi_cap_stop();
  pcu::Finalize();
  exit(exit_code);
}

void printUsage(char *argv0) {
  printf("USAGE: %s [-agwvm] <size-field> <create_file.cre>\n", argv0);
  printf("Flags:\n"
  "-a\tEvaluate size-field analytically.\n"
  "-g\tForce mesh generation.\n"
  "-v\tEnable verbose output.\n"
  "-w\tWrite before.vtk, after.vtk, and after.cre.\n"
  "-m\tConvert mesh to MDS during adaptation.\n"
  "SIZE-FIELDS:\n"
  "%d, for uniform anisotropic size-field\n"
  "%d, for wing-shock size-field\n"
  "%d, for cube-shock size-field\n"
  "%d, for cylinder boundary-layer size-field\n", 1, 2, 3, 4);
}

} // namespace

int main(int argc, char** argv) {
  // Initialize parallelism.
  pcu::Init(&argc, &argv);
  {
  pcu::PCU PCUObj;

  // Initialize logging.
  lion_set_verbosity(1);
  lion_set_stdout(stdout);
  lion_set_stderr(stderr);

  // Check arguments or print usage.
  if (argc < 3) {
    if (PCUObj.Self() == 0) {
      printUsage(argv[0]);
    }
    myExit(EXIT_FAILURE);
  }

  // Parse arguments.
  bool volume_flag = false, write_flag = false, analytic_flag = false,
       verbose_flag = false, mds_flag = false;
  for (int i = 1; i < argc - 2; ++i) {
    if (*argv[i] == '-') {
      for (int j = 1; argv[i][j] != '\0'; ++j) {
        switch(argv[i][j]) {
        case 'a':
          analytic_flag = true;
          break;
        case 'g':
          volume_flag = true;
          break;
        case 'v':
          verbose_flag = true;
          lion_set_verbosity(1);
          break;
        case 'w':
          write_flag = true;
          break;
        case 'm':
          mds_flag = true;
          break;
        default:
          printf("Error: invalid flag.\n");
          printUsage(argv[0]);
          myExit(EXIT_FAILURE);
        }
      }
    }
  }

  const char* createFileName = argv[argc - 1];
  int mode = atoi(argv[argc - 2]);

  // Initialize GMI.
  gmi_cap_start();
  gmi_register_cap();
  // create an instance of the Capstone Module activating SMLIB/CREATE/CREATE
  // for the Geometry/Mesh/Attribution databases
  const std::string gdbName("Geometry Database : SMLIB");
  const std::string mdbName("Mesh Database : Create");
  const std::string adbName("Attribution Database : Create");

  CapstoneModule  cs("capVol", gdbName.c_str(), mdbName.c_str(), adbName.c_str());

  GeometryDatabaseInterface     *g = cs.get_geometry();
  MeshDatabaseInterface         *m = cs.get_mesh();
  AppContext                    *c = cs.get_context();

  PCU_ALWAYS_ASSERT(g);
  PCU_ALWAYS_ASSERT(m);
  PCU_ALWAYS_ASSERT(c);

  // Load Capstone mesh.
  v_string filenames;
  filenames.push_back(createFileName);
  M_GModel gmodel = cs.load_files(filenames);

  M_MModel mmodel;
  if (volume_flag) {
    mmodel = cs.generate_mesh();
    if (mmodel.is_invalid()) {
      lion_eprint(1, "FATAL: Failed to mesh the model.\n");
      myExit(EXIT_FAILURE);
    }
    MG_API_CALL(m, set_current_model(mmodel));
  } else {
    // Use the first existing mesh model.
    std::vector<M_MModel> mmodels;
    MG_API_CALL(m, get_associated_mesh_models(gmodel, mmodels));
    PCU_ALWAYS_ASSERT(mmodels.size() == 1);
    mmodel = mmodels[0];
    MG_API_CALL(m, set_current_model(mmodel));
  }

  if (write_flag) {
    cs.save_file("core_capVol_before.cre", gmodel);
  }

  // Calculate adjacencies.
  MG_API_CALL(m, set_adjacency_state(REGION2FACE|REGION2EDGE|REGION2VERTEX|
				     FACE2EDGE|FACE2VERTEX));
  MG_API_CALL(m, set_reverse_states());
  MG_API_CALL(m, compute_adjacency());

  // Make APF adapter over Capstone mesh.
  ma::Mesh* apfCapMesh = apf::createMesh(m, g, &PCUObj);
  apf::writeVtkFiles("core_capVol_cap.vtk", apfCapMesh);

  ma::Mesh* adaptMesh = apfCapMesh;
  if (mds_flag) {
    adaptMesh = apf::createMdsMesh(apfCapMesh->getModel(), apfCapMesh);
    apf::reorderMdsMesh(adaptMesh);
    // APF default routine will typically fail to verify surface meshes.
    if (volume_flag) adaptMesh->verify();
    if (write_flag) {
      apf::writeVtkFiles("core_capVol_mds.vtk", adaptMesh);
    }
    apfCapMesh->setModel(nullptr); // Disown the model.
    delete apfCapMesh;
    apfCapMesh = nullptr;
  }

  // Choose appropriate size-field.
  ma::AnisotropicFunction* sf = nullptr;
  switch (mode) {
    case 1:
      sf = new UniformAniso(adaptMesh);
      break;
    case 2:
      sf = new WingShock(adaptMesh, 50);
      break;
    case 3:
      sf = new Shock(adaptMesh);
      break;
    case 4:
      sf = new CylBoundaryLayer(adaptMesh);
      break;
    default:
      lion_eprint(1, "FATAL: Invalid size-field.\n");
      myExit(EXIT_FAILURE);
  }

  // Make pumi fields for the frames and scales for anisotropic size-fields.
  apf::Field* frameField = nullptr;
  apf::Field* scaleField = nullptr;
  ma::Input *in = nullptr;
  if (!analytic_flag || write_flag) {
    frameField = apf::createFieldOn(adaptMesh, "adapt_frames", apf::MATRIX);
    scaleField = apf::createFieldOn(adaptMesh, "adapt_scales", apf::VECTOR);

    ma::Entity *v;
    ma::Iterator* it = adaptMesh->begin(0);
    while( (v = adaptMesh->iterate(it)) ) {
      ma::Vector s;
      ma::Matrix f;
      sf->getValue(v, f, s);
      apf::setVector(scaleField, v, 0, s);
      apf::setMatrix(frameField, v, 0, f);
    }
    adaptMesh->end(it);

    if (write_flag) {
      apf::writeVtkFiles("core_capVol_before.vtk", adaptMesh);
    }
  }

  if (!analytic_flag) {
    // Pass the field data.
    in = ma::makeAdvanced(ma::configure(adaptMesh, scaleField, frameField));
  } else {
    // Pass the function.
    in = ma::makeAdvanced(ma::configure(adaptMesh, sf));
  }

  ma::adapt(in);

  if (volume_flag) {
    // We can't verify surface meshes.
    adaptMesh->verify();
  }

  if (write_flag) {
    if (mds_flag) {
      apf::writeVtkFiles("core_capVol_after_mds.vtk", adaptMesh);
      MG_API_CALL(m, create_associated_model(mmodel, gmodel, "MeshAdapt"));
      MG_API_CALL(m, set_adjacency_state(REGION2FACE|REGION2EDGE|REGION2VERTEX|
                 FACE2EDGE|FACE2VERTEX));
      MG_API_CALL(m, set_reverse_states());
      // MG_API_CALL(m, compute_adjacency()); // unnecessary because no elements?
      ma::Mesh* newCapMesh = apf::createMesh(m, g);
      apf::convert(adaptMesh, newCapMesh);
      apf::writeVtkFiles("core_capVol_after_cap.vtk", newCapMesh);
      apf::destroyMesh(newCapMesh);
    } else {
      apf::writeVtkFiles("core_capVol_after_cap.vtk", adaptMesh);
    }
    cs.save_file("core_capVol_after.cre", gmodel);
  }

  /* PRINT ADAPTED MESH INFO */
  if (verbose_flag) {
    M_MModel mmdl;
    m->get_current_model(mmdl);
    std::string info;
    m->print_info(mmdl, info);
    lion_oprint(1, "%s", info.c_str());
  }

  // Clean up.
  apf::destroyMesh(adaptMesh);
  delete sf;

  // Exit calls.
  gmi_cap_stop();
  }
  pcu::Finalize();
}
