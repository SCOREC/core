#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfConvert.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <apf.h>
#include <lionPrint.h>
#include <pcu_util.h>

#include "MeshSim.h"
#include "MeshSimAdapt.h"
#include "SimDiscrete.h"
#include "SimMessages.h"
#include <iostream>

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <memory>


using namespace std;

typedef vector<double> vec;
typedef vector<vec>    mat;


void printModelStats(pGModel model);
void makeSimxModelAndMesh(
    double* coords, int nverts,
    apf::Gid*    conn,   int nelem,
    pMesh& mesh, pDiscreteModel& model,
    pVertex* vReturn, pEntity* eReturn);
bool checkVertexOrder(
    apf::Mesh2* mesh,
    pVertex* vReturn,
    int nverts);
bool checkVertexOrder(
    apf::Mesh2* mesh,
    const apf::GlobalToVert & map,
    double* coords);
void getSizesAndFrames(
    apf::Mesh2* m,
    apf::Field* sizes,
    apf::Field* frames,
    vector<apf::Vector3>& sz,
    vector<apf::Matrix3x3>& fr);
pMSAdapt addSizesToSimxMesh(
    pMesh mesh,
    int nverts,
    pVertex* vReturn,
    const vector<apf::Vector3>& sizes,
    const vector<apf::Matrix3x3>& frames);
double runSimxAdapt(pMSAdapt adapter);
double runSimxMeshImprover(
    pMesh mesh,
    double minQuality);
void destructSimxMesh(
    pMesh simxMesh,
    double*& adaptedCoords,
    int*& adaptedConns,
    int& nverts, int& nelem,
    vector<apf::Vector3> &  adaptedSizes,
    vector<apf::Matrix3x3> & adaptedFrames);
void addFields(apf::Mesh2* m,
    const apf::GlobalToVert & map,
    const char* sizeName,
    const char* frameName,
    const vector<apf::Vector3> & adaptedSizes,
    const vector<apf::Matrix3x3> &adaptedFrames);
apf::Mesh2* convertToPumi(
    pMesh simxMesh, int dim,
    const char* sizeName, const char* frameName, pcu::PCU *PCUObj);
int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  auto PCUObj = std::unique_ptr<pcu::PCU>(new pcu::PCU(MPI_COMM_WORLD));
  lion_set_verbosity(1);
  MS_init(); // Call before calling Sim_readLicenseFile
  Sim_readLicenseFile(0);
  SimDiscrete_start(0);  // initialize GeomSim Discrete library

  if (argc < 7) {
    if (PCUObj.get()->Self() == 0) {
      printf("USAGE: %s <model.dmg> <mesh.smb> <prefix>"
      	  "<scale field name> <frame field name> <min_quality>\n", argv[0]);
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  gmi_register_mesh();
  gmi_register_null();

  PCU_ALWAYS_ASSERT_VERBOSE(PCUObj.get()->Peers() == 1,
      "This utility only works for serial meshes!");

  const char* inputPumiModel = argv[1];
  const char* inputPumiMesh = argv[2];
  const char* prefix = argv[3];
  const char* sizeName = argv[4];
  const char* frameName = argv[5];
  double minQuality = atof(argv[6]);

  char outSimxModel[256];
  char outInitialSimxMesh[256];
  char outAdaptedSimxMesh[256];
  char outAdaptedPumiMesh[256];
  char outImprovedSimxMesh[256];
  char outImprovedPumiMesh[256];

  snprintf(outSimxModel, 256, "%s.smd", prefix);
  snprintf(outInitialSimxMesh, 256, "%s_initial.sms", prefix);
  snprintf(outAdaptedSimxMesh, 256, "%s_adapted.sms", prefix);
  snprintf(outAdaptedPumiMesh, 256, "%s_adapted.smb", prefix);
  snprintf(outImprovedSimxMesh, 256, "%s_adapted_improved.sms", prefix);
  snprintf(outImprovedPumiMesh, 256, "%s_adapted_improved.smb", prefix);

  apf::Mesh2* m = apf::loadMdsMesh(inputPumiModel, inputPumiMesh, PCUObj.get());

  char message[512];
  // first find the sizes field
  apf::Field* sizes  = m->findField(sizeName);
  snprintf(message, 512, "Couldn't find a field with name %s in mesh!", sizeName);
  PCU_ALWAYS_ASSERT_VERBOSE(sizes, message);

  // then find the frames field if they exist
  apf::Field* frames;
  frames  = m->findField(frameName);
  snprintf(message, 512, "Couldn't find a field with name %s in mesh!", frameName);
  PCU_ALWAYS_ASSERT_VERBOSE(frames, message);

  // remove every field except for sizes and frames
  int index = 0;
  while (m->countFields() > 2) {
    apf::Field* f = m->getField(index);
    if (f == sizes || f == frames) {
      index++;
      continue;
    }
    m->removeField(f);
    apf::destroyField(f);
  }

  m->verify();

  // extract the coordinates and connectivities
  apf::Gid* conn;
  double* coords;
  int nelem;
  int etype;
  int nverts;

  int dim = m->getDimension();
  extractCoords(m, coords, nverts);
  destruct(m, conn, nelem, etype);

  // make/save Simx mesh and model
  printf("\n===CONVERTING THE PUMI MESH TO SIMX===\n");
  pMesh simxMesh = 0;
  pDiscreteModel simxModel = 0;
  pVertex* vReturn = new pVertex[nverts];
  pEntity* eReturn = new pEntity[nelem];
  makeSimxModelAndMesh(coords, nverts, conn, nelem, simxMesh, simxModel, vReturn, eReturn);
  printf("===DONE===\n");

  printf("\n===WRITING THE SIMX MODEL %s===\n", outSimxModel);
  GM_write(simxModel,outSimxModel,0,0); // save the discrete model
  cout<<"Model stats: "<<endl;
  printModelStats(simxModel);
  printf("===DONE===\n");

  printf("\n===WRITING THE SIMX INITIAL MESH %s===\n", outInitialSimxMesh);
  M_write(simxMesh,outInitialSimxMesh, 0,0);  // write out the initial mesh data
  printf("===DONE===\n");

  printf("\n===RUNNING SIMX ADAPT===\n");
  PCU_ALWAYS_ASSERT_VERBOSE(checkVertexOrder(m, vReturn, nverts),
      "The verts orders in the pumi mesh and the created simx mesh appear to be different!\n");
  vector<apf::Vector3> sz;
  vector<apf::Matrix3x3> fr;
  getSizesAndFrames(m, sizes, frames, sz, fr);

  pMSAdapt adapter = addSizesToSimxMesh(simxMesh, nverts, vReturn, sz, fr);
  double adaptTime = runSimxAdapt(adapter);
  printf("\nSIMX ADAPT RUN-TIME: %f \n", adaptTime);
  printf("===DONE===\n");


  printf("\n===WRITING THE SIMX/SMB ADAPTED MESHES===\n");
  printf("%s\n", outAdaptedSimxMesh);
  printf("%s\n", outAdaptedPumiMesh);
  M_write(simxMesh,outAdaptedSimxMesh, 0,0);  // write out the initial mesh data
  apf::Mesh2* m2 = convertToPumi(simxMesh, dim, sizeName, frameName, PCUObj.get());

  m2->writeNative(outAdaptedPumiMesh);
  printf("===DONE===\n");

  printf("\n===RUNNING SIMX IMPROVER WITH TARGET QUALITY %f===\n", minQuality);
  double improveTime = runSimxMeshImprover(simxMesh, minQuality);
  printf("\nSIMX IMPROVER RUN-TIME: %f \n", improveTime);
  printf("===DONE===\n");

  printf("\n===WRITING THE SIMX/SMB IMPROVED MESHES===\n");
  printf("%s\n", outImprovedSimxMesh);
  printf("%s\n", outImprovedPumiMesh);
  M_write(simxMesh,outImprovedSimxMesh, 0,0);  // write out the initial mesh data
  apf::Mesh2* m3 = convertToPumi(simxMesh, dim, sizeName, frameName, PCUObj.get());

  m3->writeNative(outImprovedPumiMesh);
  printf("===DONE===\n");

    // cleanup
  M_release(simxMesh);
  GM_release(simxModel);

  m->destroyNative();
  apf::destroyMesh(m);

  m2->destroyNative();
  apf::destroyMesh(m2);

  m3->destroyNative();
  apf::destroyMesh(m3);

  SimDiscrete_stop(0);
  Sim_unregisterAllKeys();
  MS_exit();
  }
  MPI_Finalize();
}

void printModelStats(pGModel model)
{
  cout<<"  Number of model vertices: "<<GM_numVertices(model)<<endl;
  cout<<"    Vertex tags: ";
  GVIter modelVertices = GM_vertexIter(model);  // initialize the iterator
  pGVertex modelVertex;
  while( (modelVertex=GVIter_next(modelVertices)) ){ // get next vertex
    cout<<" "<<GEN_tag(modelVertex);
  }
  GVIter_delete(modelVertices);
  cout<<endl;
  cout<<"  Number of model edges: "<<GM_numEdges(model)<<endl;
  cout<<"  Number of model faces: "<<GM_numFaces(model)<<endl;
  cout<<"    Face tags: ";
  GFIter modelFaces = GM_faceIter(model);  // initialize the iterator
  pGFace modelFace;
  while( (modelFace=GFIter_next(modelFaces)) ){ // get next face
    cout<<" " << GEN_tag(modelFace);
  }
  GFIter_delete(modelFaces);
  cout<<""<<endl;
  cout<<"  Number of model regions: "<<GM_numRegions(model)<<endl;
}

void makeSimxModelAndMesh(
    double* coords, int nverts,
    apf::Gid*    connGid,   int nelem,
    pMesh& mesh, pDiscreteModel& model,
    pVertex* vReturn, pEntity* eReturn)
{
  // assuming all tet elements
  int* elementType = new int[nelem];
  for (int i = 0; i < nelem; i++)
    elementType[i] = 10;

  Sim_setMessageHandler(0);

  const int connSize = 4*nelem;
  int* conn = new int[connSize];
  for(int i=0; i<connSize; i++) {
    conn[i] = static_cast<int>(connGid[i]);
  }

  mesh = M_new(0,0);
  if(M_importFromData(mesh,nverts,coords,nelem,
      elementType,conn,vReturn,eReturn,0)) { //check for error 
    cerr<<"Error importing mesh data"<<endl;
    M_release(mesh);
    return;
  }

  // check the input mesh for intersections
  // this call must occur before the discrete model is created
  if(MS_checkMeshIntersections(mesh,0,0)) {
    cerr<<"There are intersections in the input mesh"<<endl;
    M_release(mesh);
    return;
  }

  // create the Discrete model
  model = DM_createFromMesh(mesh, 0, 0);
  if(!model) { //check for error
    cerr<<"Error creating Discrete model from mesh"<<endl;
    M_release(mesh);
    return;
  }
  DM_findEdgesByFaceNormals(model, 0, 0);
  DM_eliminateDanglingEdges(model, 0);
  if(DM_completeTopology(model, 0)) { //check for error
    cerr<<"Error completing Discrete model topology"<<endl;
    M_release(mesh);
    GM_release(model);
    return;
  }
  delete [] conn;
}

bool checkVertexOrder(
    apf::Mesh2* mesh,
    pVertex* vReturn,
    int nverts)
{
  PCU_ALWAYS_ASSERT((size_t)nverts == mesh->count(0));
  double tol = 1.e-12;
  apf::MeshEntity* ent;
  apf::MeshIterator* it;

  it = mesh->begin(0);
  int count = 0;
  while( (ent = mesh->iterate(it)) ){
    // get the coordinates for pumi vert
    apf::Vector3 pumiVertCoords;
    mesh->getPoint(ent, 0, pumiVertCoords);
    // get the coordinates for simx vert
    double xyz[3];
    pVertex vertex = vReturn[count];
    V_coord(vertex,xyz);
    apf::Vector3 simxVertCoords(xyz[0], xyz[1], xyz[2]);
    // compare the length of the difference between the two
    if ( (pumiVertCoords - simxVertCoords).getLength() > tol)
      return false;
    count++;
  }
  mesh->end(it);
  return true;
}


bool checkVertexOrder(
    apf::Mesh2* mesh,
    const apf::GlobalToVert & map,
    double* coords)
{
  double tol = 1.e-12;

  APF_CONST_ITERATE(apf::GlobalToVert, map, it)
  {
    int key = it->first;
    apf::MeshEntity* vert = it->second;
    apf::Vector3 pt;
    mesh->getPoint(vert, 0, pt);
    apf::Vector3 ptFromInput(coords[3*key],
			     coords[3*key+1],
			     coords[3*key+2]);
    if ( (pt - ptFromInput).getLength() > tol)
      return false;
  }
  return true;
}

void getSizesAndFrames(
    apf::Mesh2* m,
    apf::Field* sizes,
    apf::Field* frames,
    vector<apf::Vector3>& sz,
    vector<apf::Matrix3x3>& fr)
{
  sz.clear();
  fr.clear();

  apf::MeshEntity* ent;
  apf::MeshIterator* it;

  it = m->begin(0);
  while( (ent = m->iterate(it)) ){
    // get the coordinates for pumi vert
    apf::Vector3 sizeVector;
    apf::Matrix3x3 frameTensor;
    apf::getVector(sizes, ent, 0, sizeVector);
    apf::getMatrix(frames, ent, 0, frameTensor);
    sz.push_back(sizeVector);
    fr.push_back(frameTensor);
  }
  m->end(it);
}

pMSAdapt addSizesToSimxMesh(
    pMesh mesh,
    int nverts,
    pVertex* vReturn,
    const vector<apf::Vector3>& sizes,
    const vector<apf::Matrix3x3>& frames)
{
  PCU_ALWAYS_ASSERT_VERBOSE((size_t)nverts == sizes.size(),
      "Expecting the size of the vector to be the same as nverts!\n");

  PCU_ALWAYS_ASSERT_VERBOSE((size_t)nverts == frames.size(),
      "Expecting the size of the vector to be the same as nverts!\n");

  pMSAdapt adapter = MSA_new(mesh, 1);

  pVertex vertex;
  for (int i = 0; i < nverts; i++) {
    vertex = vReturn[i];
    apf::Vector3 sz = sizes[i];
    apf::Matrix3x3 fr = frames[i];

    apf::Vector3 row0 = apf::Vector3(fr[0][0],
				     fr[1][0],
				     fr[2][0]);
    apf::Vector3 row1 = apf::Vector3(fr[0][1],
				     fr[1][1],
				     fr[2][1]);
    apf::Vector3 row2 = apf::Vector3(fr[0][2],
				     fr[1][2],
				     fr[2][2]);

    row0 = row0 * sz[0];
    row1 = row1 * sz[1];
    row2 = row2 * sz[2];

    double anisoSize[3][3] = { {row0[0], row0[1], row0[2]},
			       {row1[0], row1[1], row1[2]},
			       {row2[0], row2[1], row2[2]} };

    MSA_setAnisoVertexSize(adapter, vertex, anisoSize);
  }

  return adapter;
}


double runSimxAdapt(pMSAdapt adapter)
{
  double t0 = pcu::Time();
  MSA_adapt(adapter, 0);
  MSA_delete(adapter);
  double t1 = pcu::Time();
  return t1 - t0;
}

double runSimxMeshImprover(pMesh mesh, double minQuality)
{
  double t0 = pcu::Time();
  pVolumeMeshImprover vmi = VolumeMeshImprover_new(mesh);
  VolumeMeshImprover_setShapeMetric(vmi, ShapeMetricType_VolLenRatio, minQuality);
  VolumeMeshImprover_execute(vmi, 0);
  VolumeMeshImprover_delete(vmi);
  double t1 = pcu::Time();
  return t1 - t0;
}

static int countVerts(pMesh mesh)
{
  VIter vit;
  pVertex vert;
  vit = M_vertexIter(mesh);
  int count = 0;
  while( (vert = VIter_next(vit)) )
    count++;
  VIter_delete(vit);
  return count;
}

static int countRegions(pMesh mesh)
{
  RIter rit;
  pRegion region;
  rit = M_regionIter(mesh);
  int count = 0;
  while( (region = RIter_next(rit)) )
    count++;
  RIter_delete(rit);
  return count;
}

static void getSizeAndFramesFromArray(
    double anisoSize[3][3],
    apf::Vector3 &sizes,
    apf::Matrix3x3 &frame)
{
  apf::Vector3 v0 = apf::Vector3(anisoSize[0][0],
				 anisoSize[0][1],
				 anisoSize[0][2]);
  apf::Vector3 v1 = apf::Vector3(anisoSize[1][0],
				 anisoSize[1][1],
				 anisoSize[1][2]);
  apf::Vector3 v2 = apf::Vector3(anisoSize[2][0],
				 anisoSize[2][1],
				 anisoSize[2][2]);
  sizes = apf::Vector3(v0.getLength(),
		       v1.getLength(),
		       v2.getLength());
  v0 = v0 / v0.getLength();
  v1 = v1 / v1.getLength();
  v2 = v2 / v2.getLength();
  /* // in pumi we store the directions in columns */
  frame = apf::Matrix3x3(v0[0], v1[0], v2[0],
			 v0[1], v1[1], v2[1],
			 v0[2], v1[2], v2[2]);
}

void destructSimxMesh(
    pMesh mesh,
    double*& adaptedCoords,
    apf::Gid*& adaptedConns,
    int& nverts, int& nelem,
    vector<apf::Vector3>&  adaptedSizes,
    vector<apf::Matrix3x3>& adaptedFrames)
{
  nverts = countVerts(mesh);
  nelem  = countRegions(mesh);

  adaptedCoords = new double[3*nverts];
  adaptedConns  = new apf::Gid[4*nelem];

  VIter vertices;
  RIter regions;
  pVertex vertex;
  pRegion region;
  double xyz[3];
  pPList regionVerts;
  int i,j;
  vertices = M_vertexIter(mesh);
  i=0;
  while( (vertex=VIter_next(vertices)) ){
    V_coord(vertex,xyz);
    adaptedCoords[i*3] = xyz[0];
    adaptedCoords[i*3+1] = xyz[1];
    adaptedCoords[i*3+2] = xyz[2];
    double anisoSize[3][3];
    V_size(vertex, NULL, anisoSize);
    apf::Vector3 sizes;
    apf::Matrix3x3 frame;
    getSizeAndFramesFromArray(anisoSize, sizes, frame);
    adaptedSizes.push_back(sizes);
    adaptedFrames.push_back(frame);
    EN_setID((pEntity)vertex,i);
    i++;
  }
  VIter_delete(vertices);
  regions = M_regionIter(mesh);
  i=0;
  while( (region = RIter_next(regions)) ){
    regionVerts = R_vertices(region,0);
    vector<apf::Gid> ids;
    for(j=0; j < 4; j++){
      vertex = (pVertex)PList_item(regionVerts,j);
      apf::Gid id = EN_id((pEntity)vertex);
      ids.push_back(id);
    }
    // simmetrix's local ordering is different from pumi's
    adaptedConns[i*4]   = ids[0];
    adaptedConns[i*4+1] = ids[2];
    adaptedConns[i*4+2] = ids[1];
    adaptedConns[i*4+3] = ids[3];
    PList_delete(regionVerts);
    i++;
  }
  RIter_delete(regions);
}

void addFields(apf::Mesh2* m,
    const apf::GlobalToVert & map,
    const char* sizeName,
    const char* frameName,
    const vector<apf::Vector3> & adaptedSizes,
    const vector<apf::Matrix3x3> &adaptedFrames)
{
  apf::Field* adaptedSizeField  = apf::createFieldOn(m, sizeName, apf::VECTOR);
  apf::Field* adaptedFrameField = apf::createFieldOn(m, frameName, apf::MATRIX);
  APF_CONST_ITERATE(apf::GlobalToVert, map, it)
  {
    int key = it->first;
    apf::MeshEntity* vert = it->second;
    apf::setVector(adaptedSizeField, vert, 0, adaptedSizes[key]);
    apf::setMatrix(adaptedFrameField, vert, 0, adaptedFrames[key]);
  }
}

apf::Mesh2* convertToPumi(
    pMesh simxMesh, int dim,
    const char* sizeName,
    const char* frameName,
    pcu::PCU *PCUObj)
{
  double* adaptedCoords;
  apf::Gid* adaptedConns;
  int adaptedNumVerts, adaptedNumElems;
  vector<apf::Vector3> adaptedSizes;
  vector<apf::Matrix3x3> adaptedFrames;
  destructSimxMesh(simxMesh, adaptedCoords, adaptedConns,
      adaptedNumVerts, adaptedNumElems, adaptedSizes, adaptedFrames);

  gmi_model* nullModel = gmi_load(".null");
  apf::Mesh2* m2 = apf::makeEmptyMdsMesh(nullModel, dim, false, PCUObj);
  apf::GlobalToVert outMap;
  apf::construct(m2, adaptedConns, adaptedNumElems, apf::Mesh::TET, outMap);;
  apf::alignMdsRemotes(m2);
  apf::deriveMdsModel(m2);
  apf::setCoords(m2, adaptedCoords, adaptedNumVerts, outMap);
  m2->verify();

  // make sure ordering of verts is what it is supposed to be.
  PCU_ALWAYS_ASSERT_VERBOSE(checkVertexOrder(m2, outMap, adaptedCoords),
      "The verts order in the adapted simx mesh and the created pumi mesh appear to be different!\n");
  // add the fields
  addFields(m2, outMap, sizeName, frameName, adaptedSizes, adaptedFrames);
  outMap.clear();
  return m2;
}
