#include <PCU.h>
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

size_t countElements(apf::Mesh2* m) {
  int dim = m->getDimension();
  size_t numElem = 0;
  for (int i = 0; i <= dim; ++i) {
    numElem += m->count(i);
  }
  return numElem;
}

int main(int argc, char** argv)
{
  /* INIT CALLS */
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();

  const char* createFileName = argv[1];

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

  std::cout << "Element count: " << countElements(apfCapMesh) << std::endl;
  std::cout << "Vertex count: " << apfCapMesh->count(0) << std::endl;

  /* SETUP AND CALL ADAPT */
  ma::Input* in;

  // choose a size field here
  ma::AnisotropicFunction* sf = new UniformAniso(apfCapMesh);

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

  /* PRINT ADAPTED MESH INFO */
  M_MModel mmdl;
  m->get_current_model(mmdl);
  std::string info;
  m->print_info(mmdl, info);
  std::cout << "Element count: " << countElements(apfCapMesh) << std::endl;
  std::cout << "Vertex count: " << apfCapMesh->count(0) << std::endl;

  long numElem = countElements(apfCapMesh), numElemRef = 5276;
  PCU_ALWAYS_ASSERT(std::abs(numElem - numElemRef) < 0.05*numElemRef);

  long numVert = apfCapMesh->count(0), numVertRef = 881;
  PCU_ALWAYS_ASSERT(std::abs(numVert - numVertRef) < 0.05*numVertRef);

  /* CLEAN UPS */
  delete sf;
  apf::destroyMesh(apfCapMesh);

  /* EXIT CALLS */
  gmi_cap_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}
