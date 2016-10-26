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
//

#include "SimPartitionedMesh.h"
#include "SimModel.h"
#include "SimUtil.h"
#include <stdlib.h>
#include <stdio.h>

#include <iostream>
#include <cstring>

using namespace std;

// Here is a simple example that loads in a serial or partitioned mesh,
// partitions it onto the supplied desired total (over all processes) number
// of partitions, and writes out the partitioned mesh. 
// To run, type mpirun -np <#procs> <exec-with-path> in parallel

void messageHandler(int type, const char *msg);

int main(int argc, char **argv)
{
  pParMesh pmesh;
  pPartitionOpts pOpts;

  // Initialize PartitionedMesh - this should be the first Simmetrix call
  // Also initializes MPI in parallel
  SimPartitionedMesh_start(&argc, &argv);

  // Read in command line arguments
  if(argc != 5){
    cout<<"Usage: "<<argv[0]<<" [model] [mesh_directory] [input_mesh] [num_parts]" << endl; // MeshDirectoryPath ModelFilename InputMeshFilename NumParts"<<endl;
    return 1;
  }

  int desiredTotNumParts = atoi(argv[4]);  // Desired total no. of partitions

  // Initialize all strings needed for filenames
  char * modelFilename = (char*) malloc((strlen(argv[1]) + 10) * sizeof(char));
  char * meshFilename = (char*) malloc((strlen(argv[2]) + strlen(argv[3]) + 10)*sizeof(char));
  char * pmeshFilename = (char*) malloc((strlen(argv[2]) + strlen(argv[3]) + 20)*sizeof(char));
  char * pmeshname = (char*) malloc((strlen(argv[3]) + 20)*sizeof(char));

  strcpy(modelFilename,argv[1]);
  strcpy( meshFilename,argv[2]);
  strcpy(pmeshFilename,argv[2]);

  //strcat(modelFilename,"/../");
  strcat( meshFilename,"/");
  strcat(pmeshFilename,"/");

  // Create filename for partitioned mesh from other input info
  memcpy(pmeshname,argv[3],strlen(argv[3])-4);
  pmeshname[strlen(argv[3])-4] = '\0';

  //strcat(modelFilename,argv[1]);
  strcat( meshFilename,argv[3]);

  strcat(pmeshFilename,pmeshname);
  strcat(pmeshFilename,"_");
  strcat(pmeshFilename,argv[4]);
  strcat(pmeshFilename,"_parts.sms");

  cout<<endl;
  cout<<"Using model and mesh:"<<endl;
  cout<<modelFilename<<endl;
  cout<< meshFilename<<endl;
  cout<<endl;
  cout<<"Partitioning into "<<argv[4]<<" parts."<<endl;  
  cout<<endl;
  cout<<"Simmetrix says..."<<endl;
  cout<<"**********************************"<<endl;


  // NOTE: Sim_readLicenseFile() is for internal testing only.  To use,
  // pass in the location of a file containing your keys.  For a release 
  // product, use Sim_registerKey() 
  Sim_readLicenseFile("/net/common/meshSim/license/license.txt");

  Sim_logOn("partition.log");

  Sim_setMessageHandler(messageHandler);
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  pGModel model = GM_load(modelFilename, 0, progress);

  // Read a serial/partitioned mesh. It's advisable to read it on a GeomSim
  // model (ie. pass valid pGModel instead of 0 below).
  pmesh = PM_load(meshFilename, sthreadNone, model, progress);

  pOpts = PM_newPartitionOpts();
  PartitionOpts_setTotalNumParts(pOpts, desiredTotNumParts);
  PM_partition(pmesh, pOpts, sthreadNone, progress);   // Do a default partitioning
  PartitionOpts_delete(pOpts);

  PM_write(pmesh, pmeshFilename, sthreadNone, progress); // Write it out to a directory
  M_release(pmesh);                                    // Delete the partitioned mesh

  cout<<"**********************************"<<endl;
  cout<<"Partitioned mesh output to:"<<endl;
  cout<<pmeshFilename<<endl;
  cout<<endl;

  Progress_delete(progress); 
  Sim_logOff();
  Sim_unregisterAllKeys();
  SimPartitionedMesh_stop();

  free(modelFilename);
  free( meshFilename);
  free(pmeshFilename);
  free(pmeshname);

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
