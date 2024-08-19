#include <PCU_C.h>
#include <apf.h>
#include <apfCAP.h>
#include <gmi.h>
#include <gmi_cap.h>
#include <apfConvert.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <ma.h>
#include <pcu_util.h>
#include <cstdlib>
#include <cstring>
#include <iostream>


#include "CapstoneModule.h"

using namespace CreateMG;
using namespace CreateMG::Attribution;
using namespace CreateMG::Mesh;
using namespace CreateMG::Geometry;

#include "capStoneSizeFields.h"


void writeCre(CapstoneModule& cs, const std::string& vtkFileName)
{
  GeometryDatabaseInterface    *gdbi = cs.get_geometry();
  MeshDatabaseInterface        *mdbi = cs.get_mesh();
  AppContext		       *ctx = cs.get_context();

  // Get the VTK writer.
  Writer *creWriter = get_writer(ctx, "Create Native Writer");
  if (!creWriter)
	  error(HERE, ERR_GENERIC, "Could not find the CRE writer");

  IdMapping idmapping;
  std::vector<M_MModel> mmodels;
  M_GModel gmodel;
  M_MModel mmodel;
  gdbi->get_current_model(gmodel);
  mdbi->get_current_model(mmodel);
  mmodels.clear();
  mmodels.push_back(mmodel);
  creWriter->write(ctx, gmodel, mmodels, vtkFileName.c_str(), idmapping);
}

int main(int argc, char** argv)
{
  /* INIT CALLS */
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();

  if (argc < 3) {
    if (PCU_Comm_Self() == 0) {
      printf("USAGE: %s <create_file.cre> <size-field>\n", argv[0]);
      printf("Size-fields:\n");
      printf("%d, for uniform anisotropic size-field\n", 1);
      printf("%d, for wing-shock size-field\n", 2);
      printf("%d, for cube-shock size-field\n", 3);
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  const char* createFileName = argv[1];
  int mode = atoi(argv[2]);

  gmi_cap_start();
  gmi_register_cap();

  /* LOAD CAPSTONE MESH */
  // create an instance of the Capstone Module activating CREATE/CREATE/CREATE
  // for the Geometry/Mesh/Attribution databases
  const std::string gdbName("Geometry Database : SMLIB");// Switch Create with SMLIB for CAD
  const std::string mdbName("Mesh Database : Create");
  const std::string adbName("Attribution Database : Create");

  CapstoneModule  cs("test", gdbName.c_str(), mdbName.c_str(), adbName.c_str());

  GeometryDatabaseInterface     *g = cs.get_geometry();
  MeshDatabaseInterface         *m = cs.get_mesh();
  AppContext                    *c = cs.get_context();

  PCU_ALWAYS_ASSERT(g);
  PCU_ALWAYS_ASSERT(m);
  PCU_ALWAYS_ASSERT(c);

  v_string filenames;
  filenames.push_back(createFileName);

  M_GModel gmodel = cs.load_files(filenames);

  M_MModel mmodel;
  // Pick the volume mesh model from associated mesh models to this geom model
  std::vector<M_MModel> mmodels;
  MG_API_CALL(m, get_associated_mesh_models(gmodel, mmodels));
  PCU_ALWAYS_ASSERT(mmodels.size() == 1);
  MG_API_CALL(m, set_current_model(mmodels[0]));

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

  /* CONVERT THE MESH TO APF::MESH2 */
  apf::Mesh2* apfCapMesh = apf::createMesh(m, g);

  /* ADD A TEST FIELD TO THE MESH TO DEMONSTRATE SOLUTION TRANSFER */
  apf::Field* tf  = apf::createFieldOn(apfCapMesh, "test_field", apf::VECTOR);
  apf::MeshEntity* ent;
  apf::MeshIterator* it = apfCapMesh->begin(0);
  while( (ent = apfCapMesh->iterate(it)) ) {
    apf::Vector3 p = ma::getPosition(apfCapMesh, ent);
    double x = p[0];
    double y = p[1];
    double z = p[2];
    apf::Vector3 s(y, z, 2*x);
    apf::setVector(tf, ent, 0, s);
  }
  apfCapMesh->end(it);

  /* WRITE THE BEFORE ADAPT MESH TO VTK USING APF VTK WRITER */
  apf::writeVtkFiles("before", apfCapMesh);


  /* SETUP AND CALL ADAPT */
  ma::Input* in;

  // choose a size field here
  ma::AnisotropicFunction* sf = 0;
  switch (mode) {
    case 1:
      sf = new UniformAniso(apfCapMesh);
      break;
    case 2:
      sf = new WingShock(apfCapMesh, 50);
      break;
    case 3:
      sf = new Shock(apfCapMesh);
      break;
  }

  // make pumi fields that hold the "frames" and "scales" for anisotropic size fields
  // here we are using user-defined size-fields. Usually, "frames" and "scales" come
  // from a solution driven error estimation procedure
  apf::Field* frameField = apf::createFieldOn(apfCapMesh, "adapt_frames", apf::MATRIX);
  apf::Field* scaleField = apf::createFieldOn(apfCapMesh, "adapt_scales", apf::VECTOR);

  apf::MeshEntity* v;
  it = apfCapMesh->begin(0);
  while( (v = apfCapMesh->iterate(it)) ) {
    apf::Vector3 s;
    apf::Matrix3x3 f;
    sf->getValue(v, f, s);
    apf::setVector(scaleField, v, 0, s);
    apf::setMatrix(frameField, v, 0, f);
  }
  apfCapMesh->end(it);

  in = ma::makeAdvanced(ma::configure(apfCapMesh, scaleField, frameField));
  in->shouldSnap = true;
  in->shouldTransferParametric = true;
  in->shouldFixShape = true;
  in->shouldForceAdaptation = true;
  if (apfCapMesh->getDimension() == 2)
    in->goodQuality = 0.04; // this is mean-ratio squared
  else // 3D meshes
    in->goodQuality = 0.027; // this is the mean-ratio cubed
  in->maximumIterations = 10;

  ma::adaptVerbose(in, false);

  /* apfCapMesh->verify(); */

  /* WRITE THE AFTER ADAPT MESH TO VTK USING APF VTK WRITER */
  apf::writeVtkFiles("after", apfCapMesh);


  /* PRINT ADAPTED MESH INFO */
  M_MModel mmdl;
  m->get_current_model(mmdl);
  std::string info;
  m->print_info(mmdl, info);

  /* WRITE THE ADAPTED MESH IN NATIVE CREATE FORMAT */
  writeCre(cs, "after.cre");

  /* CLEAN UPS */
  delete sf;

  /* EXIT CALLS */
  gmi_cap_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}
