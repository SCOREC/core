#include <PCU.h>
#include <pcu_util.h>
#include <lionPrint.h>
#include "phOutput.h"
#include "phGrowthCurves.h"
#include "phLinks.h"
#include "phAdjacent.h"
#include "apfSIM.h"
#include "gmi_sim.h"
#include <SimUtil.h>
#include <SimPartitionedMesh.h>
#include <SimAdvMeshing.h>
#include <stdlib.h>

namespace ph {
void getThinSectionStack(Output& o)
{
  o.nThinSectionStacks = 0;
  o.nThinSectionStackMeshVertices = 0;

  Input& in = *o.in;
  if (in.simmetrixMesh == 1) {
    if (in.writeSimLog)
      Sim_logOn("getThinStack.log");
    pProgress progress = Progress_new();
    Progress_setDefaultCallback(progress);

    // get simmetrix mesh
    apf::MeshSIM* apf_msim = dynamic_cast<apf::MeshSIM*>(o.mesh);
    pParMesh parMesh = apf_msim->getMesh();
    pMesh mesh = PM_mesh(parMesh,0);

    // get simmetrix model
    gmi_model* gmiModel = apf_msim->getModel();
    pGModel model = gmi_export_sim(gmiModel);

//  Algorithm: Get growth curve info
 
    typedef std::pair <pGEntity, pGFace> gPair_t;
    typedef std::multimap <pGEntity, pGFace> gPairMap_t;
    typedef std::pair <gPairMap_t::iterator, gPairMap_t::iterator> gPairMap_equalRange_t;

//  Create an empty list (gEntities) for storing gEntity
//  Create an empty multimap (gPairMap) for storing pairs gPair {KEY: gEntity, CONTENT: gFace}
//  //gEntity is the model entity where a base mesh vertex is classified
//  //gFace is the model face where 3D thin section/extrusion attribute is placed
    pPList gEntities = PList_new();
    gPairMap_t gPairs;
    gPairMap_t::iterator gPairIter;
    gPair_t gPair;

    pGEntity gEntity;
    pGFace gFace;
    pGEdge gEdge;
    pGVertex gVertex;
    pVertex vertex;
	pFace face;

    pPList gEdges = PList_new();
    pPList gVertices = PList_new();
	pPList gSrcFaces = PList_new();
	pPList regions_stack = PList_new();
	pPList faces_stack = PList_new();
	pPList gDestFaces =PList_new();

//  //generate gEntities and gPairs
//  //gEntities contains non-duplicated items
//  //gPairs may contain duplicated items
    PList_clear(gEntities);
    gPairIter = gPairs.begin();
//  FOR each model face (gFace)
    GFIter gFIter = GM_faceIter(model);
    while((gFace = GFIter_next(gFIter))){
	 if(PList_contains(gDestFaces, gFace)==0){				//check if gFace is already in gSrcFaces
//    IF gFace has extrusion/thin stack mesh attribute
  	  bool isThinStackFace = false;
      FIter fIter = M_classifiedFaceIter(mesh, gFace, 1);
      while((face = FIter_next(fIter))){
        if(EN_isExtrusionEntity(face) == 1){
		  PList_appUnique(gSrcFaces, gFace);
		  for(int j = 0; j < 1; j++){
			if((EN_isExtrusionEntity(F_region(face, j))) == 1){
				region = F_region(face, j);
				Extrusion_3DRegionsAndLayerFaces(region, regions_stack, faces_stack, 0));
				PList_appUnique(gDestFaces, F_whatIn(PList_item(faces_stack(PList_size(faces_stack)))));
				PList_clear(regions_stack);
				PList_clear(faces_stack);
            break;
			}
		  }
	    }
	  }
      GFIter_delete(gFIter);
	 }
	}
	
//  get seeds of all growth curves
    pPList allSeeds = PList_new();
    pPList gFaces = PList_new();

    pPList seeds = PList_new();
	pPList base_vertices = PList_new();
	pPList fvertices = PList_new();
	pPList allStackVertices = PList_new();
	pPList StackVertices = PList_new();
	pPList iTSnv_list

	pEntity seed;
    pGRegion gRegion;
	pRegion region;

//  FOR each gEntity in gEntities
    for(int i = 0; i < PList_size(gSrcFaces); i++){
      gEntity = (pGEntity)PList_item(gSrcFaces,i);

//    Generate a non-duplicated list (gFaces) for storing model faces associated with the key gEntity in gPairs
      PList_clear(gFaces);

//    Get mesh faces classified on gEntity
      FIter fIter = M_classifiedFaceIter(mesh, gEntity, 1);

//    FOR each face
      while((face = FIter_next(fIter))){
//      Create an empty list (seeds) for storing potential seed edges of vertex
        PList_clear(seeds);
		PList_clear(regions_stack);
		PList_clear(faces_stack);

		for(int j = 0; j < 1; j++){
			if((EN_isExtrusionEntity(F_region(face, j))) == 1){
				region = F_region(face, j);
				PList_appUnique(seeds,region); 						//add region to seeds List
				
				//Get region stack and faces_stack
				if(!(Extrusion_3DRegionsAndLayerFaces(region, regions_stack, faces_stack, 0)) == 1){
					lion_oprint(1,"%s: unexpected Extrusion_3DRegionsAndLayerFaces return value\n",__func__);
					exit(EXIT_FAILURE);
				}
				
				PList_clear(fvertices);
				
				fvertices = F_vertices((pFace)PList_item(faces_stack,0),1);					//get mesh vertices for base face
				PList_appPListUnique(base_vertices, PList(fvertices));					//add to base vertices unique list base_vertices
				
				for(int k = 0; k < PList_size(fvertices); k++){								//Loop over vertices of base face
					if(PList_contains(base_vertices, PList_item(fvertices(k)))==0){ 		//Check if vertex is already in list. If not
						PList_appUnique(base_vertices, PList_item(fvertices(k)))		//add vertex to list
						PList_clear(StackVertices);
						Plist_append(iTSnv_list, PList_size(faces_stack);
						for(int faces_stack_iter = 0; faces_stack_iter < PList_size(faces_stack); k++){					//Loop over faces in faces_stack
							PList_clear(fvertices2);
							fvertices2 = F_vertices((pFace)PList_item(faces_stack,faces_stack_iter),1);   //Get Stack vertices
							PList_append(StackVertices, PList_item(fvertices2,k));		//For face kk in connectivity stack, add vertex k to allStackVertices
						}
						PList_appPList(allStackVertices, StackVertices);
					}
				}	
				
			}
		}
	}

	int nTS = PList_size(base_vertices);
	o.nThinSectionStacks = nTS;
	
	int nTSMeshVertices = PList_size(allStackVertices)						
	o.nThinSectionStackMeshVertices = nTSMeshVertices;

	o.arrays.iTSnv = new int[nTS];

	for(int iTSnv_iter = 0; iTSnv_iter < PList_size(iTSnv_list); iTSnv_iter++){
		 o.arrays.iTSnv[i] = int(PList_item(iTSnv_list,iTSnv_iter));
	}

    o.arrays.iTSlv = new apf::MeshEntity*[nTSMeshVertices];

    for(int i = 0; i < PList_size(allStackVertices); i++){
      vertex = (pVertex)PList_item(allStackVertices,i);

      apf::MeshEntity* me = reinterpret_cast<apf::MeshEntity*> (vertex);
      o.arrays.iTSlv[i] = me;
    }

    lion_oprint(1,"%s: rank %d, ngc, nv: %d, %d\n", __func__, PCU_Comm_Self(), nTS, nTSMeshVertices);

    PCU_Add_Ints(&nTS,sizeof(nTS));
    PCU_Add_Ints(&nTSMeshVertices,sizeof(nTSMeshVertices));

    if(PCU_Comm_Self() == 0)
      lion_oprint(1,"%s: total ngc, nv: %d, %d\n", __func__, nTS, nTSMeshVertices);

	//Delete all lists
	PList_delete(gEdges);
	PList_delete(gEntities);
	PList_delete(gVertices);
	PList_delete(gEdges);
	PList_delete(allSeeds);
	PList_delete(gFaces);
	PList_delete(seeds);
	PList_delete(regions_stack);
	PList_delete(faces_stack);
	PList_delete(base_vertices);
	PList_delete(fvertices);
	PList_delete(allStackVertices);
	PList_delete(StackVertices);
	PList_delete(iTSnv_list);
	
    //clean up utility
    Progress_delete(progress);
    if (in.writeSimLog)
      Sim_logOff();
  }
  else {
    if(PCU_Comm_Self() == 0)
      lion_oprint(1,"%s: warning! not implemented for MDS mesh\n",__func__);
  }
  return;
}
} //end namespace ph