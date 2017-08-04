#include "SimModel.h"
#include "SimUtil.h"
#include "SimParasolidKrnl.h"
#include "SimPartitionedMesh.h"
#include "SimAdvMeshing.h"

#include <pcu_util.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <iostream>

using namespace std;

const char* modelFilename;
const char* meshFilename;
const char* attribFilename;
int expectedNum;

void messageHandler(int type, const char *msg);

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  SimPartitionedMesh_start(&argc, &argv);
  // Read in command line arguments
  if (argc == 5) {
    modelFilename = argv[1];
    attribFilename = argv[2];
    meshFilename = argv[3];
    expectedNum = atoi(argv[4]);
  }
  else if (argc == 4) {
    modelFilename = argv[1];
    attribFilename = argv[2];
    meshFilename = argv[3];
    expectedNum = 0;
  }
  else {
    cout<<"Usage1: "<<argv[0]<<" [model.xmt_txt] [attribute.smd] [input_mesh.sms]" << endl;
    cout<<"Usage2: "<<argv[0]<<" [model.xmt_txt] [attribute.smd] [input_mesh.sms] [expectedNum]" << endl;
    return 1;
  }

  /* print message */
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if (rank==0) {
    cout<<endl;
    cout<<"Using model and mesh: "<<modelFilename<<" "<<meshFilename<<endl;
    cout<<"Simmetrix says..."<<endl;
    cout<<"**********************************"<<endl;
  }

  Sim_readLicenseFile("/net/common/meshSim/license/license.txt");
  Sim_logOn("countBL.log");
  SimParasolid_start(1);
  SimAdvMeshing_start();

  Sim_setMessageHandler(messageHandler);
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  pNativeModel nmodel = ParasolidNM_createFromFile(modelFilename, 0);
  pGModel model = GM_load(attribFilename, nmodel, progress);
  /* ---------------------------add test code begin------------------------ */
  pParMesh pmesh = PM_load(meshFilename, model, progress);
  pMesh mesh = PM_mesh(pmesh, 0);

  // declaration
  pVertex meshVertex;
  VIter vIter;

  int counter = 0;
  vIter = M_vertexIter(mesh);
  while((meshVertex = VIter_next(vIter))) {
    if(EN_isBLEntity(meshVertex))
      counter++;
  }
  VIter_delete(vIter);

  printf("rank = %d; there are %d BL entities in this part\n", rank, counter);

  int totalNum;
  MPI_Reduce(&counter, &totalNum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if (!rank)
    PCU_ALWAYS_ASSERT(totalNum >= expectedNum);

  if (!rank)
    printf("there are %d BL entities in total\n", totalNum);

  /* --------------------------- add test code end ------------------------ */
  M_release(pmesh);
  GM_release(model);
  NM_release(nmodel);

  if (rank==0) {
    cout<<"**********************************"<<endl;
    cout<<"Done!"<<endl;
  }

  Progress_delete(progress);
  Sim_logOff();
  SimAdvMeshing_stop();
  SimParasolid_stop(1);
  Sim_unregisterAllKeys();
  SimPartitionedMesh_stop();
  MPI_Finalize();
  return 0;
}

void messageHandler(int type, const char *msg)
{
  switch (type) {
  case Sim_InfoMsg:
    fprintf(stderr,"Info: %s\n",msg);
    break;
  case Sim_DebugMsg:
    fprintf(stderr,"Debug: %s\n",msg);
    break;
  case Sim_WarningMsg:
    fprintf(stderr,"Warning: %s\n",msg);
    break;
  case Sim_ErrorMsg:
    fprintf(stderr,"Error: %s\n",msg);
    break;
  }
  return;
}

