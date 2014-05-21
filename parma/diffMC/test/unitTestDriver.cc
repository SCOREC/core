#include "parma_diffmc.h"
#include "parma_partinfo.h"
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

using apf::Mesh;
using apf::createMesh;

/**
 * @brief test parma part weight collection functions
 *
 * @return zero on success, non-zero otherwise
 */
int weightTest3D(Mesh* mesh) { 
   int priority[4] = {0, 1, 3, 2};
   priorityList pl;
   pl.sort(priority);
   MeshTag* wtag = mesh->createDoubleTag("weight",1);
   double maxImbW[4] = {60, 270, 370, 166}; // not correct
   partInfo p(mesh, wtag, pl, 0, maxImbW);
   mesh->destroyTag(wtag);
   //for(size_t i=0; i<p.size(); i++) p[i].print();
   double w[4][4] = {
      64.000,277.000,379.000,165.000,
      57.000,262.000,371.000,165.000,
      60.000,270.000,376.000,165.000,
      65.000,279.000,380.000,165.000 };
   if ( 3 != p.adjPartIds.size() ) return 1;
   for(size_t j=0; j<p.adjPartIds.size(); j++) 
      if( !isEqlArr(p.adjPartWeights[j], w[ p.adjPartIds[j] ] ) ) 
	 return 1;  
   if( !isEqlArr(p.weight, w[PCU_Comm_Self()]) ) return 1;
   return 0;
}

/**
 * @brief test parma priority sort functions
 *
 * @return zero on success, non-zero otherwise
 */
int sortTest() {
   int ierr = 0;
   {
      int priority[4] = {9, 23, 1, 6};
      priorityList pl;
      pl.sort(priority);
      if ( pl.entDim[0] != 1 || pl.priority[0] != 23 ||
	    pl.entDim[1] != 0 || pl.priority[1] != 9  ||
	    pl.entDim[2] != 3 || pl.priority[2] != 6  ||
	    pl.entDim[3] != 2 || pl.priority[3] != 1    ) {
	 ierr += 1;
      }
   }

   {
      int priority[4] = {2,1,3,1};
      priorityList pl;
      pl.sort(priority);
      if ( pl.entDim[0] != 2 || pl.priority[0] != 3 ||
	    pl.entDim[1] != 0 || pl.priority[1] != 2  ||
	    pl.entDim[2] != 3 || pl.priority[2] != 1  ||
	    pl.entDim[3] != 1 || pl.priority[3] != 1    ) {
	 ierr += 1;
      }
   }
   return ierr;
}

int main(int argc, char** argv) {
   MPI_Init(&argc,&argv);
   PCU_Comm_Init();

   int ierr = 0;
   ierr += sortTest();

  PUMI_Init(MPI_COMM_WORLD);
  pGeomMdl model = new pumi::NullModel("x");

   int commSize; 
   PCU_Comm_Size(&commSize);
   if ( 4 == commSize ) {
      pMeshMdl mesh;
      PUMI_Mesh_Create(model,mesh);
      PUMI_Mesh_LoadFromFile(mesh, "part.sms", 1);
      Mesh* mesh2 = (Mesh*) createMesh(mesh);
      ierr += weightTest3D(mesh2);
      delete mesh2;
   } 

   if ( ierr ) fprintf(stderr, "STATUS: tests failed\n");
   else fprintf(stderr, "STATUS: tests passed\n");
   MPI_Finalize();
   return ierr;
}
