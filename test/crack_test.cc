#include <apf.h>
#ifdef HAVE_SIMMETRIX
#include <SimUtil.h>
#include <MeshSim.h>
#include <gmi_sim.h>
#include <SimModel.h>
#endif
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <ma.h>
#include <crv.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <pcu_util.h>
#include <lionPrint.h>
#include <vector>
#include <cassert>
#include <stdlib.h>
#include <sstream>
#include <fstream>


void bCurver(const char* modelFile, const char* meshFile,
    pcu::PCU *PCUObj, int order, const char* outputPrefix,
    int blendOrder = 0)
{
  apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile,PCUObj);
  m->verify();
  if (m->getPCU()->Self() == 0)
    printf("attempting to curve the mesh to %d order Bezier\n", order);
  crv::BezierCurver bc(m, order, blendOrder);
  bc.run();
  m->verify();
  if (m->getPCU()->Self() == 0)
    printf("succeeded!\n");

  crv::writeCurvedVtuFiles(m, apf::Mesh::TRIANGLE, order + 2, outputPrefix);
  crv::writeCurvedWireFrame(m, order + 10, outputPrefix);

  m->destroyNative();
  apf::destroyMesh(m);
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);

  if (argc < 2) {
    if (PCUObj.Self() == 0) {
      printf("USAGE: %s <model.x_t>\n", argv[0]);
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  const char* modelFile   = argv[1];
  
  lion_set_verbosity(1);
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_null();
  gmi_register_mesh();
  gmi_model* g = gmi_load(modelFile);

  apf::Mesh2* mesh = apf::makeEmptyMdsMesh(g, 2, false, &PCUObj);

  std::vector<apf::MeshEntity*> verts;
  // vertex coordinates
  apf::Vector3 p0(-0.45, 0.3, 0.0);
  apf::Vector3 p1( 0.00, 0.3, 0.0);
  apf::Vector3 p2( 0.45, 0.3, 0.0);
  apf::Vector3 p3( 0.45,-0.3, 0.0);
  apf::Vector3 p4( 0.00,-0.3, 0.0);
  apf::Vector3 p5(-0.45,-0.3, 0.0);
  apf::Vector3 p6(-0.15, 0.0, 0.0);
  apf::Vector3 p7( 0.15, 0.0, 0.0);

  apf::ModelEntity* c;
  apf::MeshEntity* vert;
  apf::Vector3 param;

  /* --- MAKING VERTICIES --- */
  /* ------------------------ */

  // V0
  vert = 0;
  c = 0;
  c = mesh->findModelEntity(0, 18);
  PCU_ALWAYS_ASSERT(c);
  vert = mesh->createVert(c);
  PCU_ALWAYS_ASSERT(vert);
  mesh->setPoint(vert, 0, p0);
  verts.push_back(vert);

  // V1
  // this is a vertex on a model edge so the parametric
  // coordinate on the model edge must be provided.
  vert = 0;
  c = 0;
  c = mesh->findModelEntity(1, 24);
  PCU_ALWAYS_ASSERT(c);
  param = apf::Vector3(0.45, 0., 0.);
  vert = mesh->createVertex(c, p1, param);
  PCU_ALWAYS_ASSERT(vert);
  verts.push_back(vert);

  // V2
  vert = 0;
  c = 0;
  c = mesh->findModelEntity(0, 2);
  PCU_ALWAYS_ASSERT(c);
  vert = mesh->createVert(c);
  PCU_ALWAYS_ASSERT(vert);
  mesh->setPoint(vert, 0, p2);
  verts.push_back(vert);

  // V3
  vert = 0;
  c = 0;
  c = mesh->findModelEntity(0, 5);
  PCU_ALWAYS_ASSERT(c);
  vert = mesh->createVert(c);
  PCU_ALWAYS_ASSERT(vert);
  mesh->setPoint(vert, 0, p3);
  verts.push_back(vert);

  // V4
  // this is a vertex on a model edge so the parametric
  // coordinate on the model edge must be provided.
  vert = 0;
  c = 0;
  c = mesh->findModelEntity(1, 23);
  PCU_ALWAYS_ASSERT(c);
  param = apf::Vector3(0.45, 0., 0.);
  vert = mesh->createVertex(c, p4, param);
  PCU_ALWAYS_ASSERT(vert);
  verts.push_back(vert);

  // V5
  vert = 0;
  c = 0;
  c = mesh->findModelEntity(0, 16);
  PCU_ALWAYS_ASSERT(c);
  vert = mesh->createVert(c);
  PCU_ALWAYS_ASSERT(vert);
  mesh->setPoint(vert, 0, p5);
  verts.push_back(vert);

  // V6
  vert = 0;
  c = 0;
  c = mesh->findModelEntity(0, 57);
  PCU_ALWAYS_ASSERT(c);
  vert = mesh->createVert(c);
  PCU_ALWAYS_ASSERT(vert);
  mesh->setPoint(vert, 0, p6);
  verts.push_back(vert);

  // V7
  vert = 0;
  c = 0;
  c = mesh->findModelEntity(0, 54);
  PCU_ALWAYS_ASSERT(c);
  vert = mesh->createVert(c);
  PCU_ALWAYS_ASSERT(vert);
  mesh->setPoint(vert, 0, p7);
  verts.push_back(vert);


  /* --- MAKING BOUNDARY EDGES --- */
  /* ----------------------------- */

  /* Outer Boundary edges are made here.
   * The classification info must be know for each of them.
   */
  apf::MeshEntity* edgeVerts[2];
  // E0
  c = 0;
  c = mesh->findModelEntity(1, 24);
  PCU_ALWAYS_ASSERT(c);

  edgeVerts[0] = verts[1];
  edgeVerts[1] = verts[0];
  apf::buildElement(mesh, c, apf::Mesh::EDGE, edgeVerts);


  // E1
  c = 0;
  c = mesh->findModelEntity(1, 24);
  PCU_ALWAYS_ASSERT(c);
  edgeVerts[0] = verts[2];
  edgeVerts[1] = verts[1];
  apf::buildElement(mesh, c, apf::Mesh::EDGE, edgeVerts);


  // E2
  c = 0;
  c = mesh->findModelEntity(1, 6);
  PCU_ALWAYS_ASSERT(c);
  edgeVerts[0] = verts[3];
  edgeVerts[1] = verts[2];
  apf::buildElement(mesh, c, apf::Mesh::EDGE, edgeVerts);


  // E3
  c = 0;
  c = mesh->findModelEntity(1, 23);
  PCU_ALWAYS_ASSERT(c);
  edgeVerts[0] = verts[4];
  edgeVerts[1] = verts[3];
  apf::buildElement(mesh, c, apf::Mesh::EDGE, edgeVerts);

  // E4
  c = 0;
  c = mesh->findModelEntity(1, 23);
  PCU_ALWAYS_ASSERT(c);
  edgeVerts[0] = verts[5];
  edgeVerts[1] = verts[4];
  apf::buildElement(mesh, c, apf::Mesh::EDGE, edgeVerts);

  // E5
  c = 0;
  c = mesh->findModelEntity(1, 14);
  PCU_ALWAYS_ASSERT(c);
  edgeVerts[0] = verts[0];
  edgeVerts[1] = verts[5];
  apf::buildElement(mesh, c, apf::Mesh::EDGE, edgeVerts);

  /* Crack Edges are made a bit differently, they must have
   * classification info while they are made. We have to use
   * them later on when we make the faces, and hence the reason
   * for getting a pointer back when we make them.
   */
  // E6
  c = 0;
  c = mesh->findModelEntity(1, 58);
  PCU_ALWAYS_ASSERT(c);
  edgeVerts[0] = verts[7];
  edgeVerts[1] = verts[6];
  apf::MeshEntity* crackEdgeTop =
    mesh->createEntity(apf::Mesh::EDGE, c, edgeVerts);

  // E7
  c = 0;
  c = mesh->findModelEntity(1, 61);
  PCU_ALWAYS_ASSERT(c);
  edgeVerts[0] = verts[7];
  edgeVerts[1] = verts[6];
  apf::MeshEntity* crackEdgeBottom =
    mesh->createEntity(apf::Mesh::EDGE, c,  edgeVerts);

  /* --- MAKING THE FACES --- */
  /* ------------------------ */
  /* Regular Faces
   * (the ones that do not include the crack edges)
   */
  c = 0;
  c = mesh->findModelEntity(2, 26);
  PCU_ALWAYS_ASSERT(c);

  apf::MeshEntity* triVerts[3];
  // F0
  triVerts[0] = verts[0];
  triVerts[1] = verts[6];
  triVerts[2] = verts[1];
  apf::buildElement(mesh, c, apf::Mesh::TRIANGLE, triVerts);

  // F1
  triVerts[0] = verts[1];
  triVerts[1] = verts[7];
  triVerts[2] = verts[2];
  apf::buildElement(mesh, c, apf::Mesh::TRIANGLE, triVerts);


  // F2
  triVerts[0] = verts[7];
  triVerts[1] = verts[3];
  triVerts[2] = verts[2];
  apf::buildElement(mesh, c, apf::Mesh::TRIANGLE, triVerts);


  // F3
  triVerts[0] = verts[4];
  triVerts[1] = verts[3];
  triVerts[2] = verts[7];
  apf::buildElement(mesh, c, apf::Mesh::TRIANGLE, triVerts);

  // F4
  triVerts[0] = verts[5];
  triVerts[1] = verts[4];
  triVerts[2] = verts[6];
  apf::buildElement(mesh, c, apf::Mesh::TRIANGLE, triVerts);

  // F5
  triVerts[0] = verts[5];
  triVerts[1] = verts[6];
  triVerts[2] = verts[0];
  apf::buildElement(mesh, c, apf::Mesh::TRIANGLE, triVerts);

  /* Ab-Normal Faces (the ones that include the crack edges)
   * These have to be made from Edges
   */
  apf::MeshEntity* triEdges[3];
  apf::MeshEntity* ev[2];
  // F6
  triEdges[0] = crackEdgeTop;

  ev[0] = verts[7];
  ev[1] = verts[1];
  triEdges[1] = findElement(mesh, apf::Mesh::EDGE, ev);

  ev[0] = verts[1];
  ev[1] = verts[6];
  triEdges[2] = findElement(mesh, apf::Mesh::EDGE, ev);

  mesh->createEntity(apf::Mesh::TRIANGLE, c, triEdges);

  // F7
  triEdges[0] = crackEdgeBottom;

  ev[0] = verts[6];
  ev[1] = verts[4];
  triEdges[1] = findElement(mesh, apf::Mesh::EDGE, ev);

  ev[0] = verts[4];
  ev[1] = verts[7];
  triEdges[2] = findElement(mesh, apf::Mesh::EDGE, ev);

  mesh->createEntity(apf::Mesh::TRIANGLE, c, triEdges);

  mesh->acceptChanges();
  mesh->verify(); // make sure the mesh is valid!

  mesh->writeNative("crack_linear.smb");
  writeVtkFiles("crack_linear", mesh);

  mesh->destroyNative();
  apf::destroyMesh(mesh);

  /* --- CURVING THE FACES --- */
  /* ------------------------- */
  for (int i = 2; i < 7; i++) {
    char output[256];
    snprintf(output,256,"crack_curved_to_order_%d", i);
    bCurver(modelFile, "crack_linear.smb", &PCUObj, i, output);
  }

#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  }
  MPI_Finalize();
}
