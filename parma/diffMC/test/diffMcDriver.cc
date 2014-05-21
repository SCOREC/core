#include "parma_diffmc.h"
#include "parma_commons.h"

#include "apf.h"
#include "apfMesh.h"
#include "apfPUMI.h"
#include "PCU.h"
#include "pumi.h"
#include "pumi_mesh.h"
#include "pumi_geom.h"
#include "modeler.h"
#include "NullModel.h"

#include <unistd.h>

using apf::Mesh;
using apf::MeshEntity;
using apf::MeshTag;
using apf::createMesh;

#define PARMA_FAIL(message)\
{fprintf(stderr,"ParMA ERROR: %s: "message"\n",__func__);\
abort();}
#define PARMA_FAIL_IF(condition,message)\
if (condition)\
PARMA_FAIL(message)

void unitWeightMesh(Mesh* m, MeshTag* wtag) {
   double w = 1;
   MeshEntity* e;
   for(int d=0; d<=m->getDimension(); d++) {
      apf::MeshIterator* itr = m->begin(d);
      while( (e = m->iterate(itr)) ) {
	 m->setDoubleTag(e, wtag, &w);
      }
      m->end(itr);
   }
}

void weightMeshForDiffusionFromSinglePart(Mesh* m, MeshTag* wtag, const int heavyPartId = 0) {
   double w = 1;
   if ( heavyPartId == PCU_Comm_Self() ) {
      w = 2;
   } 
   MeshEntity* e;
   for(int d=0; d<=m->getDimension(); d++) {
      apf::MeshIterator* itr = m->begin(d);
      while( (e = m->iterate(itr)) ) {
	 m->setDoubleTag(e, wtag, &w);
      }
      m->end(itr);
   }
}

void weightMeshForDiffusionToSinglePart(Mesh* m, MeshTag* wtag, const int heavyPartId = 0) {
   double w = 1;
   if ( heavyPartId != PCU_Comm_Self() ) {
      w = 2;
   }
   MeshEntity* e;
   for(int d=0; d<=m->getDimension(); d++) {
      apf::MeshIterator* itr = m->begin(d);
      while( (e = m->iterate(itr)) ) {
	 m->setDoubleTag(e, wtag, &w);
      }
      m->end(itr);
   }
}

void broadcastDouble(double& val) {
    MPI_Bcast(&val, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
void broadcastInt(int& val) {
    MPI_Bcast(&val, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
void broadcastIntArray(int* arr, int len) {
    MPI_Bcast(arr, len, MPI_INT, 0, MPI_COMM_WORLD);
}

void broadcastString(char* src, const size_t maxStringLength) {
    MPI_Bcast(src, maxStringLength, MPI_CHAR, 0, MPI_COMM_WORLD);
}

typedef struct inputArgs {
   typedef int (*IntArrayFour)[4];
   int maxStringLen;
   char* meshFileName;
   int testNum;
   int (*entPriority)[4];
   int maxItr;
   double stallTol;
   int dbgLvl;
   void buffPrint(char* msg) {
    sprintf(msg, "Inputs: meshFileName=%s testNum=%d\n",
            meshFileName, testNum);
   }
   void print() {
     int rank;
     PCU_Comm_Rank(&rank);

     char msg[1024];
     buffPrint(msg);
     if ( 0 == rank) printf(msg);
   }
   void usage() {
     int rank;
     PCU_Comm_Rank(&rank);
     if ( 0 == rank) {
	fprintf(stderr, "Usage: <exe> -m meshFileName -t testNum\n"
	                "-m specify the input mesh file name\n"
	                "-t specify test number to run\n"
	                "-i specify the number of iterations per entity type\n"
	                "-[vefr] <prio> specify the priority of [v]tx, [e]dge, [f]ace, [r]egion. The larger the integer the higher the priority.\n"
	                "-s specify the stall tolerance between 0%% and 100%%\n"
	                "-d enable debug\n");
     }
   }
   inputArgs() {
      maxStringLen = 2048;
      meshFileName = new char[maxStringLen];
      strcpy(meshFileName, "");
      testNum = 0;
      entPriority = (IntArrayFour) malloc(sizeof(*entPriority));
      for (int i=0; i<4; i++) 
         (*entPriority)[i] = 0;
      maxItr = 20;
      stallTol = 0.5;
      dbgLvl = 0;
   }
   ~inputArgs() {
      free(entPriority);
      delete [] meshFileName;
   }
   void read(int argc, char** argv) {
      int rank;
      PCU_Comm_Rank(&rank);

      int opt;
      if (0 == rank) {
	 while ((opt = getopt(argc, argv, "m:t:v:e:f:r:i:s:d")) != -1) {
	    switch (opt) {
	       case 'm':
		  strncpy(meshFileName, optarg, maxStringLen);
		  meshFileName[maxStringLen - 1] = '\0';
		  break;
	       case 't':
		  testNum = atoi(optarg);
		  break;
               case 'v':
		  (*entPriority)[0] = atoi(optarg);
		  break;                    
	       case 'e':
		  (*entPriority)[1] = atoi(optarg);
		  break;
	       case 'f':
		  (*entPriority)[2] = atoi(optarg);
		  break;
	       case 'r':
		  (*entPriority)[3] = atoi(optarg);
		  break;              
	       case 'i':
		  maxItr = atoi(optarg);
		  break;              
	       case 's':
		  stallTol = atof(optarg);
		  break;              
	       case 'd':
		  dbgLvl = 1;
		  break;              
	       default:
                  usage();
		  SCUTIL_Finalize();
		  exit(EXIT_FAILURE);
	    }
	 }
      }
      broadcastString(meshFileName, maxStringLen);
      broadcastInt(testNum);
      broadcastIntArray(*entPriority, 4);
      broadcastInt(maxItr);
      broadcastDouble(stallTol);
      broadcastInt(dbgLvl);
   }
   void check() {
      int rank;
      PCU_Comm_Rank(&rank);
      if ( 0 == strcmp(meshFileName, "") && 0 == testNum ) {
	 if (rank == 0) usage();
	 PARMA_FAIL("input mesh file name and testNum required");
      }
   }
} inArgs;

void printTestUsage() {
   fprintf(stderr, 
           "t=0 user inputs\n"
           "t=1 size(part==0) = 2*size(part!=0), improve rgn balance\n"
           "t=2 size(part==0) = 2*size(part!=0), improve vtx balance\n"
           "t=3 size(part==0) = 2*size(part!=0), improve rgn>vtx balance\n"
           "t=4 size(part!=0) = 2*size(part==0), improve rgn balance\n");
}

void runTest(Mesh* m, inArgs& args) {
   switch(args.testNum) {
      case 0: {
                 MeshTag* wtag = m->createDoubleTag("weight",1);
                 unitWeightMesh(m, wtag);
		 Parma parma(m, wtag);
		 int ierr = parma.run(args.entPriority, args.dbgLvl, args.maxItr, 
                               args.stallTol, 1.05);
		 break;
	      }
      case 1: {
                 MeshTag* wtag = m->createDoubleTag("weight",1);
                 weightMeshForDiffusionFromSinglePart(m, wtag);
		 int priority[4] = {0,0,0,1};
		 Parma parma(m, wtag);
		 int ierr = parma.run(&priority, args.dbgLvl, args.maxItr, 
                               args.stallTol, 1.05);
                 PARMA_FAIL_IF(0 != ierr, "diffusion test failed");
		 break;
	      }
      case 2: {
                 MeshTag* wtag = m->createDoubleTag("weight",1);
                 weightMeshForDiffusionFromSinglePart(m, wtag);
		 int priority[4] = {1,0,0,0};
		 Parma parma(m, wtag);
		 int ierr = parma.run(&priority, args.dbgLvl, args.maxItr, 
                               args.stallTol, 1.05);
                 PARMA_FAIL_IF(0 != ierr, "diffusion test failed");
		 break;
	      }
      case 3: {
                 MeshTag* wtag = m->createDoubleTag("weight",1);
                 weightMeshForDiffusionFromSinglePart(m, wtag);
		 int priority[4] = {1,0,0,2};
		 Parma parma(m, wtag);
		 int ierr = parma.run(&priority, args.dbgLvl, args.maxItr, 
                               args.stallTol, 1.05);
                 PARMA_FAIL_IF(0 != ierr, "diffusion test failed");
		 break;
	      }
      case 4: {
                 MeshTag* wtag = m->createDoubleTag("weight",1);
                 weightMeshForDiffusionToSinglePart(m, wtag);
		 int priority[4] = {0,0,0,1};
		 Parma parma(m, wtag);
		 int ierr = parma.run(&priority, args.dbgLvl, args.maxItr, 
                               args.stallTol, 1.05);
                 PARMA_FAIL_IF(0 != ierr, "diffusion test failed");
		 break;
	      }
      default: {
		  if( 0 == PCU_Comm_Self() ) printTestUsage();
		  PARMA_FAIL("Invalid test number specified... exiting");
	       }
   }
}

void loadMesh(pGeomMdl& model, pMeshMdl& mesh, inArgs& args, Mesh*& m) {
   model = new pumi::NullModel("x");
   PUMI_Mesh_Create(model,mesh);
   PUMI_Mesh_LoadFromFile(mesh, args.meshFileName, 1);

   int isValid = 0;
   PUMI_Mesh_Verify(mesh,&isValid);
   PARMA_FAIL_IF(0 == isValid, "input mesh is not valid");
   

   m = createMesh(mesh);
}

int main(int argc, char** argv) {
   PUMI_Init(MPI_COMM_WORLD);
  
   inArgs args; 
   args.read(argc, argv); 
   args.print();
   args.check();

   pGeomMdl model;
   pMeshMdl pumiMesh;

   Mesh* mesh;
   loadMesh(model, pumiMesh, args, mesh);

   //PUMI_Mesh_WriteToFile(pumiMesh, "in.vtk", PCU_Comm_Peers()-1);

   runTest(mesh, args);

   //PUMI_Mesh_WriteToFile(pumiMesh, "dcFinal.vtk", PCU_Comm_Peers()-1);

   delete mesh;

   int isValid;
   PUMI_Mesh_Verify(pumiMesh,&isValid);
   //PUMI_Mesh_WriteToFile(pumiMesh, "out.vtk", PCU_Comm_Peers()-1);
   PUMI_Mesh_Del(pumiMesh);
   if (PCU_Comm_Self()) 
      PUMI_Geom_Del(model);
   PUMI_Finalize();
   return 0;
}
