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

  MS_init(); // initial MeshSim library
  SimModel_start();
  Sim_readLicenseFile("/net/common/meshSim/license/license.txt");
  Sim_logOn("createDM.log");
  SimParasolid_start(1);
  SimDiscrete_start(0);

  Sim_setMessageHandler(messageHandler);
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  pNativeModel nmodel = ParasolidNM_createFromFile(modelFilename, 0);
  pGModel model = GM_load(attribFilename, nmodel, progress);
  pMesh mesh = M_load(meshFilename, model, progress);
  pDiscreteModel dmodel = 0;
  NM_release(nmodel);

  // Print out info of mesh
  cout<<"Input mesh: "<<endl;
  cout<<"Number of mesh vertices: "<<M_numVertices(mesh)<<endl;
  cout<<"Number of mesh edges: "<<M_numEdges(mesh)<<endl;
  cout<<"Number of mesh faces: "<<M_numFaces(mesh)<<endl;
  cout<<"Number of mesh regions: "<<M_numRegions(mesh)<<endl;

  // load discrete model from parasolid model and mesh
  dmodel = DM_createFromModel(model,mesh);

  // release the original mesh
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

  // load discrete model again
  dmodel = (pDiscreteModel) GM_load(outmodelFilename, 0, progress);
  mesh = DM_getMeshCopy(dmodel,0);

  // Print out info of mesh
  cout<<"Internal mesh: "<<endl;
  cout<<"Number of mesh vertices: "<<M_numVertices(mesh)<<endl;
  cout<<"Number of mesh edges: "<<M_numEdges(mesh)<<endl;
  cout<<"Number of mesh faces: "<<M_numFaces(mesh)<<endl;
  cout<<"Number of mesh regions: "<<M_numRegions(mesh)<<endl;

  cout<<"load and check again Discrete model: "<<endl;
  cout<<"Number of model vertices: "<<GM_numVertices(dmodel)<<endl;
  cout<<"Number of model edges: "<<GM_numEdges(dmodel)<<endl;
  cout<<"Number of model faces: "<<GM_numFaces(dmodel)<<endl;
  cout<<"Number of model regions: "<<GM_numRegions(dmodel)<<endl;

  M_release(mesh);
  GM_release(dmodel);
  Progress_delete(progress);
  SimDiscrete_stop(0);
  SimParasolid_stop(1);
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
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

