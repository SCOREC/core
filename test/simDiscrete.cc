#ifdef HAVE_SIMMETRIX
#include "MeshSim.h"
#include "SimDiscrete.h"
#include "SimMessages.h"
#include "SimParasolidKrnl.h"
#endif

#include <iostream>

using namespace std;

const char* modelFilename;
const char* attribFilename;
const char* meshFilename;
const char* outmodelFilename;

void messageHandler(int type, const char *msg);

int main(int argc, char *argv[])
{
  // read in command line arguments
  if (argc == 5) {
    modelFilename    = argv[1];
    attribFilename   = argv[2];
	meshFilename     = argv[3];
	outmodelFilename = argv[4];
  }
  else {
    cout<<"Usage:"<<argv[0]<<" [model.x_t][model.smd][mesh.sms][output_model.smd]"<<endl;
	return 1;
  }

  cout<<"Using model and mesh: "<<modelFilename<<" "<<meshFilename<<endl;

  Sim_readLicenseFile("/net/common/meshSim/license/license.txt");
  Sim_logOn("createDM.log");
  SimUtil_start();
  MS_init(); // initial MeshSim library
  SimParasolid_start(1);
  SimDiscrete_start(0);

  Sim_setMessageHandler(messageHandler);
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  pNativeModel nmodel = ParasolidNM_createFromFile(modelFilename, 0);
  pGModel model = GM_load(attribFilename, nmodel, progress);
  pMesh mesh = M_load(meshFilename, model, progress);
  pDiscreteModel dmodel = 0;

  // check the input mesh for intersections
  // this call must occur before the discrete model is created
  if(MS_checkMeshIntersections(mesh,0,progress)) {
    cerr<<"There are intersections in the input mesh"<<endl;
    M_release(mesh);
    GM_release(model);
	NM_release(nmodel);
    return 1;
  }

  // create the Discrete model
  dmodel = DM_createFromMesh(mesh, 1, progress);
  if(!dmodel) { //check for error
    cerr<<"Error creating Discrete model from mesh"<<endl;
    M_release(mesh);
    GM_release(model);
	NM_release(nmodel);
    return 1;
  }

  // define the Discrete model
  DM_findEdgesByFaceNormals(dmodel, 20, progress);
  DM_eliminateDanglingEdges(dmodel, progress);
  if(DM_completeTopology(dmodel, progress)) { //check for error
    cerr<<"Error completing Discrete model topology"<<endl;
    M_release(mesh);
    GM_release(model);
    GM_release(dmodel);
	NM_release(nmodel);
    return 1;
  }

  // Since we told the Discrete model to use the input mesh, we release our
  // pointer to it.  It will be fully released when the Discrete model is released.
  M_release(mesh);

  // Print out information about the model
  cout<<"Parametric model: "<<endl;
  cout<<"Number of model vertices: "<<GM_numVertices(model)<<endl;
  cout<<"Number of model edges: "<<GM_numEdges(model)<<endl;
  cout<<"Number of model faces: "<<GM_numFaces(model)<<endl;
  cout<<"Number of model regions: "<<GM_numRegions(model)<<endl;

  cout<<"Discrete model: "<<endl;
  cout<<"Number of model vertices: "<<GM_numVertices(dmodel)<<endl;
  cout<<"Number of model edges: "<<GM_numEdges(dmodel)<<endl;
  cout<<"Number of model faces: "<<GM_numFaces(dmodel)<<endl;
  cout<<"Number of model regions: "<<GM_numRegions(dmodel)<<endl;

  GM_write(dmodel,outmodelFilename,0,progress); // save the discrete model
  GM_release(model);
  GM_release(dmodel);
  NM_release(nmodel);
  Progress_delete(progress);
  SimDiscrete_stop(0);
  SimParasolid_stop(1);
  MS_exit();
  Sim_unregisterAllKeys();
  SimUtil_stop();
  Sim_logOff();
}

void messageHandler(int type, const char *msg)
{
  switch (type) {
  case Sim_InfoMsg:
    cout<<"Info: "<<msg<<endl;
    break;
  case Sim_DebugMsg:
    cout<<"Debug: "<<msg<<endl;
    break;
  case Sim_WarningMsg:
    cout<<"Warning: "<<msg<<endl;
    break;
  case Sim_ErrorMsg:
    cout<<"Error: "<<msg<<endl;
    break;
  }
  return;
}

