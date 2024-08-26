#include <apf.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <pcu_util.h>
#include <apfDynamicVector.h>
#include <apfDynamicMatrix.h>
#include <crv.h>
#include <cassert>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif

// === includes for safe_mkdir ===
#include <reel.h>
#include <sys/types.h> /*required for mode_t for mkdir on some systems*/
#include <sys/stat.h> /*using POSIX mkdir call for SMB "foo/" path*/
#include <errno.h> /* for checking the error from mkdir */


static void safe_mkdir(
    const char* path);

static void makeCavityMeshes(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    apf::Mesh2* &cavityMeshLinear,
    apf::Mesh2* &cavityMeshCurved);

static void makeEntMeshes(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    apf::Mesh2* &entMeshLinear,
    apf::Mesh2* &entMeshCurved);

static void writeMeshes(
    apf::Mesh2* m,
    const char* prefix0,
    const char* prefix1,
    const char* prefix2,
    const char* name,
    int res = 10);

static apf::MeshTag* tagMesh(
    apf::Mesh2* m,
    int dim,
    int model);

static apf::MeshTag* tagMesh(
    apf::Mesh2* m,
    const std::vector<std::string>& ids);

int main(int argc, char** argv)
{

  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  if (PCUObj.Peers() > 1) {
    printf("%s should only be used for serial (single part) meshes!\n", argv[0]);
    printf("use the serialize utility to get a serial mesh, and retry!\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  if (argc != 6) {
    printf("USAGE: %s <model> <mesh> <prefix> <resolution> <mode>\n", argv[0]);
    printf("modes are as follows \n");
    printf("aa: creates all vert, edge, face cavities\n");
    printf("ai: creates all vert, edge, face cavities classified on interior\n");
    printf("ab: creates all vert, edge, face cavities classified on boundary\n");
    printf("va: creates all vert cavities\n");
    printf("vi: creates all vert cavities classified on interior\n");
    printf("vb: creates all vert cavities classified on boundary\n");
    printf("ea: creates all edge cavities\n");
    printf("ei: creates all edge cavities classified on interior\n");
    printf("eb: creates all edge cavities classified on boundary\n");
    printf("fa: creates all face cavities\n");
    printf("fi: creates all face cavities classified on interior\n");
    printf("fb: creates all face cavities classified on boundary\n");
    printf("ls: get a list from user and creates cavities for that list\n");
    printf("tagname: creates cavities for all entities that have tagname\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif

  gmi_register_null();

  const char* modelFile = argv[1];
  const char* meshFile  = argv[2];
  const char* prefix    = argv[3];
  int         res       = atoi(argv[4]);
  std::string mode(argv[5]);

  apf::MeshTag* tag = 0;


  // load the mesh and check if the tag exists on the mesh
  apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile,&PCUObj);

  if (mode.compare(std::string("aa")) == 0)
    tag = tagMesh(m, -1, 1);
  else if (mode.compare(std::string("ai")) == 0)
    tag = tagMesh(m, -1, 2);
  else if (mode.compare(std::string("ab")) == 0)
    tag = tagMesh(m, -1, 3);
  else if (mode.compare(std::string("va")) == 0)
    tag = tagMesh(m, 0, 1);
  else if (mode.compare(std::string("vi")) == 0)
    tag = tagMesh(m, 0, 2);
  else if (mode.compare(std::string("vb")) == 0)
    tag = tagMesh(m, 0, 3);
  else if (mode.compare(std::string("ea")) == 0)
    tag = tagMesh(m, 1, 1);
  else if (mode.compare(std::string("ei")) == 0)
    tag = tagMesh(m, 1, 2);
  else if (mode.compare(std::string("eb")) == 0)
    tag = tagMesh(m, 1, 3);
  else if (mode.compare(std::string("fa")) == 0)
    tag = tagMesh(m, 2, 1);
  else if (mode.compare(std::string("fi")) == 0)
    tag = tagMesh(m, 2, 2);
  else if (mode.compare(std::string("fb")) == 0)
    tag = tagMesh(m, 2, 3);
  else if (mode.compare(std::string("ls")) == 0) {
    std::cout << "provide the list of entities format \"v_i,e_i,f_i\"" << std::endl;
    std::cout << "example: v120 e23 f0 represents vert 120," << std::endl;
    std::cout << "edge 23 and face 0." << std::endl;
    std::cout << "total #verts=" << m->count(0);
    std::cout << ", #edges=" << m->count(1);
    std::cout << ", #faces=" << m->count(2) << std::endl;
    int cnt=0;
    std::cout << "enter #of ents in the list:" << std::endl;
    std::cin >> cnt;
    std::vector<std::string> ents;
    ents.clear();
    while ((int)ents.size() < cnt) {
      std::string temp;
      std::cin >> temp;
      ents.push_back(temp);
    }
    PCU_ALWAYS_ASSERT((int)ents.size() == cnt);
    std::cout << "creating cavities for " << std::endl;
    for (int i = 0; i < (int)ents.size(); i++)
      std::cout << ents[i] << " ";
    std::cout << std::endl;

    tag = tagMesh(m, ents);
  }
  else {
    tag = m->findTag(mode.c_str());
    if (!tag) {
      printf("tag with name %s was not found on the mesh. Aborting!\n", mode.c_str());
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
  }

  PCU_ALWAYS_ASSERT(tag);


  // make the root directory to save the cavity info
  safe_mkdir(prefix);
  writeMeshes(m, prefix, "mesh", NULL, "linear", res);
  // change the order of the mesh
  m->changeShape(crv::getBezier(3), true);
  writeMeshes(m, prefix, "mesh", NULL, "curved", res);


  apf::MeshEntity* e;
  apf::MeshIterator* it;


  for (int d = 0; d < 3; d++) {
    it = m->begin(d);
    int index = 0;

    // for now cavities are defined as follows
    // all the upward adjacent entities of dimension "dim" that
    // are adjacent to the verts of the "e". E.g., in case of edges,
    // this would give us the bi-directional edge collapse cavity.
    while ( (e = m->iterate(it)) ) {
      if (!m->hasTag(e, tag)) {
      	index++;
      	continue;
      }
      int etype = m->getType(e);

      apf::Mesh2* cavityMeshCurved = 0;
      apf::Mesh2* cavityMeshLinear = 0;
      makeCavityMeshes(m, e, cavityMeshLinear, cavityMeshCurved);

      apf::Mesh2* entMeshCurved = 0;
      apf::Mesh2* entMeshLinear = 0;
      makeEntMeshes(m, e, entMeshLinear, entMeshCurved);

      char cavityFolderName[128];
      char cavityFileNameLinear[128];
      char cavityFileNameCurved[128];
      char entityFileNameLinear[128];
      char entityFileNameCurved[128];
      snprintf(cavityFolderName, 128, "%s_%05d", apf::Mesh::typeName[etype], index);
      snprintf(cavityFileNameLinear, 128, "%s", "cavity_linear");
      snprintf(cavityFileNameCurved, 128, "%s", "cavity_curved");
      snprintf(entityFileNameLinear, 128, "%s", "entity_linear");
      snprintf(entityFileNameCurved, 128, "%s", "entity_curved");
      writeMeshes(cavityMeshLinear, prefix, "cavities",
      	  cavityFolderName, cavityFileNameLinear, res);
      writeMeshes(cavityMeshCurved, prefix, "cavities",
      	  cavityFolderName, cavityFileNameCurved, res);

      writeMeshes(entMeshLinear, prefix, "cavities",
      	  cavityFolderName, entityFileNameLinear, res);
      writeMeshes(entMeshCurved, prefix, "cavities",
      	  cavityFolderName, entityFileNameCurved, res);

      entMeshLinear->destroyNative();
      entMeshCurved->destroyNative();
      apf::destroyMesh(entMeshLinear);
      apf::destroyMesh(entMeshCurved);
      cavityMeshLinear->destroyNative();
      cavityMeshCurved->destroyNative();
      apf::destroyMesh(cavityMeshLinear);
      apf::destroyMesh(cavityMeshCurved);
      index++;
    }
    m->end(it);
  }

  // rest of the clean up
  m->destroyNative();
  apf::destroyMesh(m);

#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif

  }
  MPI_Finalize();
}

static void safe_mkdir(
    const char* path)
{
  mode_t const mode = S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH;
  int err;
  errno = 0;
  err = mkdir(path, mode);
  if (err != 0 && errno != EEXIST)
  {
    reel_fail("Err: could not create directory \"%s\"\n", path);
  }
}

static apf::Vector3 getEdgeCenter(
    apf::Mesh2* m,
    apf::MeshEntity* e)
{
  PCU_ALWAYS_ASSERT(m->getType(e) == apf::Mesh::EDGE);
  apf::MeshEntity* dv[2];
  m->getDownward(e, 0, dv);
  apf::Vector3 center(0., 0., 0.);
  for (int i = 0; i < 2; i++) {
    apf::Vector3 p;
    m->getPoint(dv[i], 0, p);
    center = center + p;
  }
  center = center * (0.5);
  return center;
}

static apf::Mesh2*  makePoint(
    apf::Mesh2* m,
    apf::MeshEntity* e)
{
  PCU_ALWAYS_ASSERT(m->getType(e) == apf::Mesh::VERTEX);
  apf::Mesh2* sphMesh = apf::makeEmptyMdsMesh(gmi_load(".null"), 1, false, m->getPCU());
  double xrange[2] = {1.e16, -1.e16};
  double yrange[2] = {1.e16, -1.e16};
  double zrange[2] = {1.e16, -1.e16};

  apf::Adjacent adj;
  m->getAdjacent(e, 1, adj);

  for (int i = 0; i < (int)adj.getSize(); i++) {
    apf::Vector3 p = getEdgeCenter(m, adj[i]);
    if (p[0] < xrange[0]) xrange[0] = p[0];
    if (p[0] > xrange[1]) xrange[1] = p[0];

    if (p[1] < yrange[0]) yrange[0] = p[1];
    if (p[1] > yrange[1]) yrange[1] = p[1];

    if (p[2] < zrange[0]) zrange[0] = p[2];
    if (p[2] > zrange[1]) zrange[1] = p[2];
  }

  double minsize = 1.e16;
  if (xrange[1]-xrange[0] < minsize) minsize = xrange[1] - xrange[0];
  if (yrange[1]-yrange[0] < minsize) minsize = yrange[1] - yrange[0];
  if (zrange[1]-zrange[0] < minsize) minsize = zrange[1] - zrange[0];

  double radius = minsize / 20.;
  apf::Vector3 center;
  m->getPoint(e, 0, center);

  // z = 0 plane
  int n = 20;
  const double pi = 3.141595;
  const apf::Vector3 param(0., 0., 0.);
  std::vector<apf::MeshEntity*> vs;
  vs.clear();
  for (int i = 0; i < n; i++) {
    apf::Vector3 p(0., 0., 0.);
    p[0] = center[0] + radius * std::cos(2.*i*pi/n);
    p[1] = center[1] + radius * std::sin(2.*i*pi/n);
    p[2] = center[2];
    apf::MeshEntity* newV = sphMesh->createVertex(m->toModel(e), p, param);
    vs.push_back(newV);
  }

  for (int i = 0; i < n; i++) {
    apf::MeshEntity* dv[2];
    dv[0] = vs[i];
    dv[1] = vs[(i+1)%20];
    sphMesh->createEntity(apf::Mesh::EDGE, m->toModel(e), dv);
  }

  // y = 0 plane
  vs.clear();
  for (int i = 0; i < n; i++) {
    apf::Vector3 p(0., 0., 0.);
    p[0] = center[0] + radius * std::cos(2.*i*pi/n);
    p[1] = center[1];
    p[2] = center[2] + radius * std::sin(2.*i*pi/n);
    apf::MeshEntity* newV = sphMesh->createVertex(m->toModel(e), p, param);
    vs.push_back(newV);
  }

  for (int i = 0; i < n; i++) {
    apf::MeshEntity* dv[2];
    dv[0] = vs[i];
    dv[1] = vs[(i+1)%20];
    sphMesh->createEntity(apf::Mesh::EDGE, m->toModel(e), dv);
  }

  // x = 0 plane
  vs.clear();
  for (int i = 0; i < n; i++) {
    apf::Vector3 p(0., 0., 0.);
    p[0] = center[0];
    p[1] = center[1] + radius * std::cos(2.*i*pi/n);
    p[2] = center[2] + radius * std::sin(2.*i*pi/n);
    apf::MeshEntity* newV = sphMesh->createVertex(m->toModel(e), p, param);
    vs.push_back(newV);
  }

  for (int i = 0; i < n; i++) {
    apf::MeshEntity* dv[2];
    dv[0] = vs[i];
    dv[1] = vs[(i+1)%20];
    sphMesh->createEntity(apf::Mesh::EDGE, m->toModel(e), dv);
  }
  sphMesh->acceptChanges();
  apf::deriveMdsModel(sphMesh);
  sphMesh->verify();
  return sphMesh;
}

static void makeEntMeshes(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    apf::Mesh2* &entMeshLinear,
    apf::Mesh2* &entMeshCurved)
{
  PCU_ALWAYS_ASSERT(!entMeshLinear);
  PCU_ALWAYS_ASSERT(!entMeshCurved);

  const apf::Vector3 param(0., 0., 0.);

  if (m->getType(e) == apf::Mesh::VERTEX)
  {
    entMeshLinear = makePoint(m, e);
    entMeshCurved = makePoint(m, e);
    return;
  }
  if (m->getType(e) == apf::Mesh::EDGE)
  {
    entMeshLinear = apf::makeEmptyMdsMesh(gmi_load(".null"), 1, false, m->getPCU());
    entMeshCurved = apf::makeEmptyMdsMesh(gmi_load(".null"), 1, false, m->getPCU());
    apf::MeshEntity* vs[2];
    m->getDownward(e, 0, vs);
    apf::Vector3 p[2];
    m->getPoint(vs[0], 0, p[0]);
    m->getPoint(vs[1], 0, p[1]);
    apf::MeshEntity* newVs[2];
    newVs[0] = entMeshLinear->createVertex(0, p[0], param);
    newVs[1] = entMeshLinear->createVertex(0, p[1], param);
    entMeshLinear->createEntity(apf::Mesh::EDGE, 0, newVs);

    apf::MeshEntity* newVsc[2];
    newVsc[0] = entMeshCurved->createVertex(0, p[0], param);
    newVsc[1] = entMeshCurved->createVertex(0, p[1], param);
    apf::MeshEntity* edge = entMeshCurved->createEntity(apf::Mesh::EDGE, 0, newVsc);

    entMeshLinear->acceptChanges();
    apf::deriveMdsModel(entMeshLinear);
    entMeshLinear->verify();

    entMeshCurved->acceptChanges();
    apf::deriveMdsModel(entMeshCurved);
    entMeshCurved->verify();

    apf::FieldShape* fs = m->getShape();
    entMeshCurved->changeShape(fs, true);
    if (fs->countNodesOn(apf::Mesh::EDGE))
    {
      for (int i = 0; i < fs->countNodesOn(apf::Mesh::EDGE); i++) {
	apf::Vector3 p;
	m->getPoint(e, i, p);
	entMeshCurved->setPoint(edge, i, p);
      }
    }
    entMeshCurved->acceptChanges();
    return;
  }
  if (m->getType(e) == apf::Mesh::TRIANGLE)
  {
    entMeshLinear = apf::makeEmptyMdsMesh(gmi_load(".null"), 2, false, m->getPCU());
    entMeshCurved = apf::makeEmptyMdsMesh(gmi_load(".null"), 2, false, m->getPCU());
    apf::MeshEntity* downverts[3];
    apf::MeshEntity* downedges[3];
    m->getDownward(e, 0, downverts);
    m->getDownward(e, 1, downedges);
    int edge_vert[3][2];
    for (int i = 0; i < 3; i++) {
      apf::MeshEntity* dv[2];
      m->getDownward(downedges[i], 0, dv);
      int i0 = apf::findIn(downverts, 3, dv[0]);
      int i1 = apf::findIn(downverts, 3, dv[1]);
      PCU_ALWAYS_ASSERT(i0 != -1);
      PCU_ALWAYS_ASSERT(i1 != -1);
      edge_vert[i][0] = i0;
      edge_vert[i][1] = i1;
    }

    apf::MeshEntity* newvertsLinear[3];
    apf::MeshEntity* newedgesLinear[3];
    apf::MeshEntity* newvertsCurved[3];
    apf::MeshEntity* newedgesCurved[3];
    apf::Vector3 param(0.,0.,0.);
    for (int i = 0; i < 3; i++) {
      apf::Vector3 p;
      m->getPoint(downverts[i], 0, p);
      newvertsLinear[i] = entMeshLinear->createVertex(0, p, param);
      newvertsCurved[i] = entMeshCurved->createVertex(0, p, param);
    }

    for (int i = 0; i < 3; i++) {
      apf::MeshEntity* evLinear[2] = {
      	newvertsLinear[edge_vert[i][0]],
      	newvertsLinear[edge_vert[i][1]]
      };
      newedgesLinear[i] = entMeshLinear->createEntity(
    	apf::Mesh::EDGE, 0, evLinear);

      apf::MeshEntity* evCurved[2] = {
      	newvertsCurved[edge_vert[i][0]],
      	newvertsCurved[edge_vert[i][1]]
      };
      newedgesCurved[i] = entMeshCurved->createEntity(
    	apf::Mesh::EDGE, 0, evCurved);
    }


    entMeshLinear->createEntity(
    	apf::Mesh::TRIANGLE, 0, newedgesLinear);

    apf::MeshEntity* face =
    entMeshCurved->createEntity(
    	apf::Mesh::TRIANGLE, 0, newedgesCurved);


    entMeshLinear->acceptChanges();
    apf::deriveMdsModel(entMeshLinear);
    entMeshLinear->verify();

    entMeshCurved->acceptChanges();
    apf::deriveMdsModel(entMeshCurved);
    entMeshCurved->verify();

    apf::FieldShape* fs = m->getShape();
    entMeshCurved->changeShape(fs, true);

    int nnodes = fs->countNodesOn(apf::Mesh::EDGE);
    if (nnodes) {
      for (int i = 0; i < 3; i++) {
	for (int n = 0; n < nnodes; n++) {
	  apf::Vector3 p;
	  m->getPoint(downedges[i], n, p);
	  entMeshCurved->setPoint(newedgesCurved[i], n, p);
	}
      }
    }

    nnodes = fs->countNodesOn(apf::Mesh::TRIANGLE);
    if (nnodes) {
      for (int n = 0; n < nnodes; n++) {
	apf::Vector3 p;
	m->getPoint(e, n, p);
	entMeshCurved->setPoint(face, n, p);
      }
    }

    entMeshCurved->acceptChanges();
    return;
  }
}

static void makeCavityMeshes(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    apf::Mesh2* &cavityMeshLinear,
    apf::Mesh2* &cavityMeshCurved)
{
  typedef std::vector<apf::MeshEntity*> Cavity;
  typedef std::vector<apf::MeshEntity*>::iterator CavityIter;
  PCU_ALWAYS_ASSERT(!cavityMeshLinear);
  PCU_ALWAYS_ASSERT(!cavityMeshCurved);

  int dim = m->getDimension();

  // input cavities
  Cavity icavity3; // tets in the cavity
  Cavity icavity2; // faces in the cavity
  Cavity icavity1; // edges in the cavity
  Cavity icavity0; // verts in the cavity
  icavity3.clear();
  icavity2.clear();
  icavity1.clear();
  icavity0.clear();


  // connectivity tables
  std::vector<std::vector<int>> tet_face;
  std::vector<std::vector<int>> face_edge;
  std::vector<std::vector<int>> edge_vert;

  apf::Downward dv;
  int nv = m->getDownward(e, 0, dv);
  for (int i = 0; i < nv; i++) {
    apf::Adjacent adj;
    m->getAdjacent(dv[i], dim, adj);
    for (int j = 0; j < (int)adj.getSize(); j++) {
      if (std::find(icavity3.begin(), icavity3.end(), adj[j]) == icavity3.end())
      	icavity3.push_back(adj[j]);
    }
  }


  // construct the unique vector of faces and tet_face connectivity table
  for (int i = 0; i < (int)icavity3.size(); i++) {
    apf::Downward dents;
    int nents = m->getDownward(icavity3[i], 2, dents);
    for (int j = 0; j < nents; j++) {
      if (std::find(icavity2.begin(), icavity2.end(), dents[j]) == icavity2.end())
      	icavity2.push_back(dents[j]);
    }
    std::vector<int> conn;
    for (int j = 0; j < nents; j++) {
      CavityIter it = std::find(icavity2.begin(), icavity2.end(), dents[j]);
      PCU_ALWAYS_ASSERT(it != icavity2.end());
      conn.push_back(std::distance(icavity2.begin(), it));
    }
    PCU_ALWAYS_ASSERT((int)conn.size() == nents);
    tet_face.push_back(conn);
  }
  PCU_ALWAYS_ASSERT(icavity3.size() == tet_face.size());

  // construct the unique vector of edges and face_edge connectivity table
  for (int i = 0; i < (int)icavity2.size(); i++) {
    apf::Downward dents;
    int nents = m->getDownward(icavity2[i], 1, dents);
    for (int j = 0; j < nents; j++) {
      if (std::find(icavity1.begin(), icavity1.end(), dents[j]) == icavity1.end())
      	icavity1.push_back(dents[j]);
    }
    std::vector<int> conn;
    for (int j = 0; j < nents; j++) {
      CavityIter it = std::find(icavity1.begin(), icavity1.end(), dents[j]);
      PCU_ALWAYS_ASSERT(it != icavity1.end());
      conn.push_back(std::distance(icavity1.begin(), it));
    }
    PCU_ALWAYS_ASSERT((int)conn.size() == nents);
    face_edge.push_back(conn);
  }
  PCU_ALWAYS_ASSERT(icavity2.size() == face_edge.size());

  // construct the unique vector of verts and edge_vert connectivity table
  for (int i = 0; i < (int)icavity1.size(); i++) {
    apf::Downward dents;
    int nents = m->getDownward(icavity1[i], 0, dents);
    for (int j = 0; j < nents; j++) {
      if (std::find(icavity0.begin(), icavity0.end(), dents[j]) == icavity0.end())
      	icavity0.push_back(dents[j]);
    }
    std::vector<int> conn;
    for (int j = 0; j < nents; j++) {
      CavityIter it = std::find(icavity0.begin(), icavity0.end(), dents[j]);
      PCU_ALWAYS_ASSERT(it != icavity0.end());
      conn.push_back(std::distance(icavity0.begin(), it));
    }
    PCU_ALWAYS_ASSERT((int)conn.size() == nents);
    edge_vert.push_back(conn);
  }
  PCU_ALWAYS_ASSERT(icavity1.size() == edge_vert.size());


  // output cavities
  Cavity ocavity3linear; // tets in the cavity
  Cavity ocavity2linear; // faces in the cavity
  Cavity ocavity1linear; // edges in the cavity
  Cavity ocavity0linear; // verts in the cavity
  Cavity ocavity3curved; // tets in the cavity
  Cavity ocavity2curved; // faces in the cavity
  Cavity ocavity1curved; // edges in the cavity
  Cavity ocavity0curved; // verts in the cavity
  ocavity3linear.clear();
  ocavity2linear.clear();
  ocavity1linear.clear();
  ocavity0linear.clear();
  ocavity3curved.clear();
  ocavity2curved.clear();
  ocavity1curved.clear();
  ocavity0curved.clear();



  // construct the cavity meshes
  // we do this in a bottom up fashion, ie vets first then edges, faces and tets
  /* cavityMeshLinear = apf::makeEmptyMdsMesh(m->getModel(), dim, false); */
  /* cavityMeshCurved = apf::makeEmptyMdsMesh(m->getModel(), dim, false); */
  cavityMeshLinear = apf::makeEmptyMdsMesh(gmi_load(".null"), dim, false, m->getPCU());
  cavityMeshCurved = apf::makeEmptyMdsMesh(gmi_load(".null"), dim, false, m->getPCU());

  // verts
  for (int i = 0; i < (int) icavity0.size(); i++) {
    apf::MeshEntity* ent = icavity0[i];
    apf::ModelEntity* c = m->toModel(ent);
    apf::Vector3 coords;
    apf::Vector3 params;
    m->getPoint(ent, 0, coords);
    m->getParam(ent, params);

    apf::MeshEntity* newEnt = cavityMeshLinear->createVertex(c, coords, params);
    ocavity0linear.push_back(newEnt);

    apf::MeshEntity* newEntc = cavityMeshCurved->createVertex(c, coords, params);
    ocavity0curved.push_back(newEntc);
  }
  PCU_ALWAYS_ASSERT(icavity0.size() == ocavity0linear.size());
  PCU_ALWAYS_ASSERT(icavity0.size() == ocavity0curved.size());

  // edges
  for (int i = 0; i < (int) icavity1.size(); i++) {
    apf::MeshEntity* ent = icavity1[i];
    apf::ModelEntity* c = m->toModel(ent);
    apf::MeshEntity* downv[2];
    downv[0] = ocavity0linear[edge_vert[i][0]];
    downv[1] = ocavity0linear[edge_vert[i][1]];
    apf::MeshEntity* newEnt = cavityMeshLinear->createEntity(
    	apf::Mesh::EDGE, c, downv);
    ocavity1linear.push_back(newEnt);

    downv[0] = ocavity0curved[edge_vert[i][0]];
    downv[1] = ocavity0curved[edge_vert[i][1]];
    apf::MeshEntity* newEntc = cavityMeshCurved->createEntity(
    	apf::Mesh::EDGE, c, downv);
    ocavity1curved.push_back(newEntc);
  }
  PCU_ALWAYS_ASSERT(icavity1.size() == ocavity1linear.size());
  PCU_ALWAYS_ASSERT(icavity1.size() == ocavity1curved.size());

  // faces
  for (int i = 0; i < (int) icavity2.size(); i++) {
    apf::MeshEntity* ent = icavity2[i];
    apf::ModelEntity* c = m->toModel(ent);
    apf::MeshEntity* downe[3];
    downe[0] = ocavity1linear[face_edge[i][0]];
    downe[1] = ocavity1linear[face_edge[i][1]];
    downe[2] = ocavity1linear[face_edge[i][2]];
    apf::MeshEntity* newEnt = cavityMeshLinear->createEntity(
    	apf::Mesh::TRIANGLE, c, downe);
    ocavity2linear.push_back(newEnt);

    downe[0] = ocavity1curved[face_edge[i][0]];
    downe[1] = ocavity1curved[face_edge[i][1]];
    downe[2] = ocavity1curved[face_edge[i][2]];
    apf::MeshEntity* newEntc = cavityMeshCurved->createEntity(
    	apf::Mesh::TRIANGLE, c, downe);
    ocavity2curved.push_back(newEntc);
  }
  PCU_ALWAYS_ASSERT(icavity2.size() == ocavity2linear.size());
  PCU_ALWAYS_ASSERT(icavity2.size() == ocavity2curved.size());


  // tets
  for (int i = 0; i < (int) icavity3.size(); i++) {
    apf::MeshEntity* ent = icavity3[i];
    apf::ModelEntity* c = m->toModel(ent);
    apf::MeshEntity* downf[4];
    downf[0] = ocavity2linear[tet_face[i][0]];
    downf[1] = ocavity2linear[tet_face[i][1]];
    downf[2] = ocavity2linear[tet_face[i][2]];
    downf[3] = ocavity2linear[tet_face[i][3]];
    apf::MeshEntity* newEnt = cavityMeshLinear->createEntity(
    	apf::Mesh::TET, c, downf);
    ocavity3linear.push_back(newEnt);

    downf[0] = ocavity2curved[tet_face[i][0]];
    downf[1] = ocavity2curved[tet_face[i][1]];
    downf[2] = ocavity2curved[tet_face[i][2]];
    downf[3] = ocavity2curved[tet_face[i][3]];
    apf::MeshEntity* newEntc = cavityMeshCurved->createEntity(
    	apf::Mesh::TET, c, downf);
    ocavity3curved.push_back(newEntc);
  }
  PCU_ALWAYS_ASSERT(icavity3.size() == ocavity3linear.size());
  PCU_ALWAYS_ASSERT(icavity3.size() == ocavity3curved.size());


  cavityMeshLinear->acceptChanges();
  apf::deriveMdsModel(cavityMeshLinear);
  cavityMeshLinear->verify();

  cavityMeshCurved->acceptChanges();
  apf::deriveMdsModel(cavityMeshCurved);
  cavityMeshCurved->verify();

  // curve cavityMeshCurved
  // this can be done by
  // a) setting the shape of cavityMeshCurved to that of m
  // b) setting the node coordinates of the entities
  // in ocavity{1,2,3}curved to those of entities
  // in icavity{1,2,3}
  cavityMeshCurved->changeShape(m->getShape(), true);
  apf::FieldShape* fs = cavityMeshCurved->getShape();

  int nnodes = fs->countNodesOn(apf::Mesh::TET);
  if (nnodes) {
    for (int i = 0; i < (int)icavity3.size(); i++) {
      apf::MeshEntity* fromEnt = icavity3[i];
      apf::MeshEntity* toEnt   = ocavity3curved[i];
      for (int j = 0; j < nnodes; j++) {
	apf::Vector3 p;
	m->getPoint(fromEnt, j, p);
	cavityMeshCurved->setPoint(toEnt, j, p);
      }
    }
  }

  nnodes = fs->countNodesOn(apf::Mesh::TRIANGLE);
  if (nnodes) {
    for (int i = 0; i < (int)icavity2.size(); i++) {
      apf::MeshEntity* fromEnt = icavity2[i];
      apf::MeshEntity* toEnt   = ocavity2curved[i];
      for (int j = 0; j < nnodes; j++) {
	apf::Vector3 p;
	m->getPoint(fromEnt, j, p);
	cavityMeshCurved->setPoint(toEnt, j, p);
      }
    }
  }

  nnodes = fs->countNodesOn(apf::Mesh::EDGE);
  if (nnodes) {
    for (int i = 0; i < (int)icavity1.size(); i++) {
      apf::MeshEntity* fromEnt = icavity1[i];
      apf::MeshEntity* toEnt   = ocavity1curved[i];
      for (int j = 0; j < nnodes; j++) {
	apf::Vector3 p;
	m->getPoint(fromEnt, j, p);
	cavityMeshCurved->setPoint(toEnt, j, p);
      }
    }
  }

  cavityMeshCurved->acceptChanges();
}

static void writeMeshes(
    apf::Mesh2* m,
    const char* prefix0,
    const char* prefix1,
    const char* prefix2,
    const char* name,
    int res)
{
  PCU_ALWAYS_ASSERT(prefix0);
  PCU_ALWAYS_ASSERT(name);
  if (!prefix1)
    PCU_ALWAYS_ASSERT(!prefix2);

  int order = m->getShape()->getOrder();
  std::stringstream ss;
  ss << prefix0 << "/";
  if (prefix1) {
    ss << prefix1;
    safe_mkdir(ss.str().c_str());
    ss << "/";
  }
  if (prefix2) {
    ss << prefix2;
    safe_mkdir(ss.str().c_str());
    ss << "/";
  }
  ss << name;

  if (order == 1) {
    apf::writeVtkFiles(ss.str().c_str(), m);
  }
  else {
    crv::writeCurvedVtuFiles(m, apf::Mesh::TRIANGLE, res, ss.str().c_str());
    crv::writeCurvedWireFrame(m, res, ss.str().c_str());
  }
  ss << ".smb";
  m->writeNative(ss.str().c_str());
}

static apf::MeshTag* tagMesh(
    apf::Mesh2* m,
    int dim,
    int model)
{
  // model = 1 means all
  // model = 2 means only interior
  // model = 3 means only boundary
  // dim = -1 tags all dims
  // dim =  0 tags verts only
  // dim =  1 tags edges only
  // dim =  2 tags faces only
  apf::MeshEntity* e;
  apf::MeshIterator* it;
  apf::MeshTag* t = m->createIntTag("which_ent", 1);
  for (int d = 0; d < 3; d++) {
    if (dim == 0 && d != 0) continue;
    if (dim == 1 && d != 1) continue;
    if (dim == 2 && d != 2) continue;
    it = m->begin(d);
    while ( (e = m->iterate(it)) ) {
      int mtype = m->getModelType(m->toModel(e));
      if (model == 2 && mtype !=3) continue;
      if (model == 3 && mtype ==3) continue;
      int val = 1; // the value does not matter
      m->setIntTag(e, t, &val);
    }
    m->end(it);
  }
  return t;
}

static void getEntIds(
    const std::vector<std::string>& ids,
    std::vector<int>& vids,
    std::vector<int>& eids,
    std::vector<int>& fids)
{
  vids.clear();
  eids.clear();
  fids.clear();

  for (std::size_t i = 0; i < ids.size(); i++) {
    std::string key = ids[i].substr(0,1);
    int value = atoi(ids[i].substr(1).c_str());
    if (key.compare(std::string("v")) == 0)
      vids.push_back(value);
    if (key.compare(std::string("e")) == 0)
      eids.push_back(value);
    if (key.compare(std::string("f")) == 0)
      fids.push_back(value);
  }

}

static apf::MeshTag* tagMesh(
    apf::Mesh2* m,
    const std::vector<std::string>& ids)
{
  std::vector<int> vids;
  std::vector<int> eids;
  std::vector<int> fids;
  getEntIds(ids, vids, eids, fids);
  PCU_ALWAYS_ASSERT(ids.size() == vids.size()+eids.size()+fids.size());

  std::vector<int>::iterator vit;
  apf::MeshEntity* e;
  apf::MeshIterator* it;
  apf::MeshTag* t = m->createIntTag("which_ent", 1);
  for (int d = 0; d < 3; d++) {
    int index = 0;
    it = m->begin(d);
    while ( (e = m->iterate(it)) ) {
      bool found = false;
      if (d == 0 &&
      	  std::find(vids.begin(), vids.end(), index) != vids.end())
      	found = true;
      if (d == 1 &&
      	  std::find(eids.begin(), eids.end(), index) != eids.end())
      	found = true;
      if (d == 2 &&
      	  std::find(fids.begin(), fids.end(), index) != fids.end())
      	found = true;
      int val = 1; // the value does not matter
      if (found)
      	m->setIntTag(e, t, &val);
      index++;
    }
    m->end(it);
  }
  return t;
}
