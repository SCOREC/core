//
// This code for making partitioned meshes using simmetrix utilities is
// modified from 
// http://www.scorec.rpi.edu/~cwsmith/SCOREC/SimModSuite_latest/PartitionedMesh/group__PMEx__Simple.html
// 
// This code takes the model, mesh and desired number of parts as input.
// Also the full path of the directory that the mesh file is in
// It assums the model file is up one directory.
// It places the partitioned mesh in the same folder as the original mesh.
// The filename of the PMesh is based on the filename of the original mesh.
//
// This can run in serial or parallel
// Running in serial corresponds to a max of 128 parts
// To get more parts, run in serial then in parallel,the filenames 
// won't work that great though so some tweaking is necessary.
// E.g. to get 256 parts, create a 2 part mesh first. Then run this
// with 'mpirun -np 2' and load the 2 part mesh and specify 256
// as the desired partitions.
//
// Dan Fovargue - Feb 2014
// Fan Yang     - Nov 2016
//

/* cheap hackish way to get SIM_PARASOLID and SIM_ACIS */
#include <PCU.h>
#include "gmi_sim_config.h"
#include <gmi_sim.h>

#ifdef HAVE_SIMMETRIX
#include "SimPartitionedMesh.h"
#include "SimModel.h"
#include "SimUtil.h"
#include "SimDiscrete.h"
#ifdef SIM_PARASOLID
#include "SimParasolidKrnl.h"
#endif
#ifdef HAVE_SIMADVMESHING
#include <SimAdvMeshing.h>
#endif
#endif

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cstring>
#include <sstream>

using namespace std;

// Here is a simple example that loads in a serial or partitioned mesh,
// partitions it onto the supplied desired total (over all processes) number
// of partitions, and writes out the partitioned mesh. 
// To run, type mpirun -np <#procs> <exec-with-path> in parallel

const char* modelFilename;
#ifdef SIM_PARASOLID
const char* attribFilename;
#endif
const char* meshFilename;
const char* outmeshFilename;
int desiredTotNumParts;  // Desired total no. of partitions

void messageHandler(int type, const char *msg);

int main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  PCU_Protect();
  // Initialize PartitionedMesh - this should be the first Simmetrix call
  // Also initializes MPI in parallel
  SimPartitionedMesh_start(&argc, &argv);
  SimDiscrete_start(0);
#ifdef HAVE_SIMADVMESHING
  SimAdvMeshing_start();
#endif
#ifdef SIM_PARASOLID
  SimParasolid_start(1);
#endif
  // Read in command line arguments
  if (argc == 4) {
    modelFilename = argv[1];
    meshFilename = argv[2];
    desiredTotNumParts = atoi(argv[3]);  // Desired total no. of partitions
  }
#ifdef SIM_PARASOLID
  else if (argc == 5) {
    modelFilename = argv[1];
    attribFilename = argv[2];
    meshFilename = argv[3];
    desiredTotNumParts = atoi(argv[4]);  // Desired total no. of partitions
  }
#endif
  else {
    cout<<"Usage1: "<<argv[0]<<" [model.smd] [input_mesh.sms] [num_parts]" << endl;
    cout<<"Enable Parasolid: " << endl;
    cout<<"Usage2: "<<argv[0]<<" [model.x_t] [attrib.smd] [input_mesh.sms] [num_parts]" << endl;
    return 1;
  }

  std::stringstream ss;
  ss << "outmesh_" << desiredTotNumParts << "_parts.sms";
  std::string tmp = ss.str();
  outmeshFilename = tmp.c_str();

  /* print message */
  if (PCU_Comm_Self()==0) {
    cout<<endl;
    cout<<"Using model and mesh: "<<modelFilename<<" "<<meshFilename<<endl;
    cout<<"Partitioning into "<< desiredTotNumParts <<" parts."<<endl;
    cout<<"Simmetrix says..."<<endl;
    cout<<"**********************************"<<endl;
  }

  // NOTE: Sim_readLicenseFile() is for internal testing only.  To use,
  // pass in the location of a file containing your keys.  For a release 
  // product, use Sim_registerKey() 
  Sim_readLicenseFile("/net/common/meshSim/license/license.txt");
  Sim_logOn("partition.log");

  Sim_setMessageHandler(messageHandler);
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  pGModel model;
#ifdef SIM_PARASOLID
  pNativeModel nmodel = NULL;
  if (argc == 5) {
    nmodel = ParasolidNM_createFromFile(modelFilename, 0);
    model = GM_load(attribFilename, nmodel, progress);
  }
  else
#endif
  {
    model = GM_load(modelFilename, 0, progress);
  }

  // Read a serial/partitioned mesh. It's advisable to read it on a GeomSim
  // model (ie. pass valid pGModel instead of 0 below).
  pParMesh pmesh = PM_load(meshFilename, model, progress);
  pPartitionOpts pOpts = PartitionOpts_new();

  PartitionOpts_setTotalNumParts(pOpts, desiredTotNumParts);
  PM_partition(pmesh, pOpts, progress);   // Do a default partitioning
  PartitionOpts_delete(pOpts);

  PM_write(pmesh, outmeshFilename, progress); // Write it out to a directory
  M_release(pmesh);                                    // Delete the partitioned mesh
  GM_release(model);
#ifdef SIM_PARASOLID
  NM_release(nmodel);
#endif

  if (PCU_Comm_Self()==0) {
    cout<<"**********************************"<<endl;
    cout<<"Partitioned mesh output to: "<<outmeshFilename<<endl;
    cout<<endl;
  }

  Progress_delete(progress);
  Sim_logOff();
#ifdef SIM_PARASOLID
  SimParasolid_stop(1);
#endif
  Sim_unregisterAllKeys();
#ifdef HAVE_SIMADVMESHING
  SimAdvMeshing_stop();
#endif
  SimDiscrete_stop(0);
  SimPartitionedMesh_stop();
  PCU_Comm_Free();
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
