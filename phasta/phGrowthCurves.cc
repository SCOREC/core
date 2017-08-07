#include <PCU.h>
#include <pcu_util.h>
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
void getGrowthCurves(Output& o)
{
  o.nGrowthCurves = 0;
  o.nLayeredMeshVertices = 0;

  Input& in = *o.in;
  if (in.simmetrixMesh == 1) {
    Sim_logOn("getGrowthCurves.log");
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
//  //gFace is the model face where 3D boundary layer attribute is placed
    pPList gEntities = PList_new();
    gPairMap_t gPairs;
    gPairMap_t::iterator gPairIter;
    gPair_t gPair;

    pGEntity gEntity;
    pGFace gFace;
    pGEdge gEdge;
    pGVertex gVertex;
    pVertex vertex;

    pPList gEdges = PList_new();
    pPList gVertices = PList_new();

//  //generate gEntities and gPairs
//  //gEntities contains non-duplicated items
//  //gPairs may contain duplicated items
    PList_clear(gEntities);
    gPairIter = gPairs.begin();
//  FOR each model face (gFace)
    GFIter gFIter = GM_faceIter(model);
    while((gFace = GFIter_next(gFIter))){
//    IF gFace has 3D boundary layer attribute
  	  bool isBoundaryLayerFace = false;
      VIter vIter = M_classifiedVertexIter(mesh, gFace, 1);
      while((vertex = VIter_next(vIter))){
        if(BL_isBaseEntity(vertex,gFace) == 1){
          isBoundaryLayerFace = true;
          break;
        }
      }

      if(isBoundaryLayerFace){
//      Add gFace to gEntities
//  		Add gPair {gFace, gFace} to gPairs
        PList_appUnique(gEntities,gFace);
        gPair = std::make_pair(gFace,gFace);
        gPairIter = gPairs.insert(gPairIter,gPair);

//      FOR each model edge (gEdge) on the closure of gFace
        gEdges = GF_edges(gFace);
        for(int i = 0; i < PList_size(gEdges); i++){
//  	    Add gEdge to gEntities
//  		  Add gPair {gEdge, gFace} to gPairs
          gEdge = (pGEdge)PList_item(gEdges,i);
          PList_appUnique(gEntities,gEdge);
          gPair = std::make_pair(gEdge,gFace);
          gPairIter = gPairs.insert(gPairIter,gPair);

//  	    FOR each model vertex (gVertex) on the closure of gEdge
          gVertices = GE_vertices(gEdge);
          for(int j = 0; j < PList_size(gVertices); j++){
//  		    Add gVertex to gEntities
//  			  Add gPair {gVertex, gFace} to gPairs
  				  gVertex = (pGVertex)PList_item(gVertices,j);
            PList_appUnique(gEntities,gVertex);
            gPair = std::make_pair(gVertex,gFace);
            gPairIter = gPairs.insert(gPairIter,gPair);
    	    }
        }
      }
    }

//  get seeds of all growth curves
    pPList allSeeds = PList_new();
    pPList gFaces = PList_new();
    gPairMap_equalRange_t gPair_equalRange;

    pPList seeds = PList_new();
    pPList blendSeeds = PList_new();
    pEntity seed;
    pGRegion gRegion;

//  FOR each gEntity in gEntities
    for(int i = 0; i < PList_size(gEntities); i++){
      gEntity = (pGEntity)PList_item(gEntities,i);

//    Generate a non-duplicated list (gFaces) for storing model faces associated with the key gEntity in gPairs
      PList_clear(gFaces);
      gPair_equalRange = gPairs.equal_range(gEntity);
      for(gPairIter = gPair_equalRange.first; gPairIter != gPair_equalRange.second; gPairIter++){
        gFace = gPairIter->second;
        PList_appUnique(gFaces,gFace);
      }

//    Get mesh vertices (vertices) classified on gEntity excluding the closure
      VIter vIter = M_classifiedVertexIter(mesh, gEntity, 0);

//    FOR each vertex in vertices
      while((vertex = VIter_next(vIter))){
//      Create an empty list (seeds) for storing potential seed edges of vertex
        PList_clear(seeds);

//      FOR each gFace in gFaces
        for(int j = 0; j < PList_size(gFaces); j++){
          gFace = (pGFace)(PList_item(gFaces,j));
//        FOR each side (faceSide) of gFace where a model region (gRegion) exists
          for(int faceSide = 0; faceSide < 2; faceSide++){
            if(!(gRegion = GF_region(gFace,faceSide)))
              continue;

            if(BL_isBaseEntity(vertex,gFace) == 0)
              continue;

            int hasSeed = BL_stackSeedEntity(vertex, gFace, faceSide, gRegion, &seed);

            switch(hasSeed){
              case 1:
                //there is 1 seed edge
                PList_appUnique(seeds,seed);
                break;
              case -1:
                //this is a blend, there will be multiple seeds
                PList_clear(blendSeeds);
                if(!(BL_blendSeedEdges(vertex, gFace, faceSide, gRegion, blendSeeds) == 1)){
                  printf("%s: unexpected BL_blendSeedEdges return value\n",__func__);
                  exit(EXIT_FAILURE);
                }
                PList_appPListUnique(seeds, blendSeeds);
                break;
              case 0:
                //there is no seed edge
                break;
              default:
                printf("%s: unexpected BL_stackSeedEntity return value\n",__func__);
                exit(EXIT_FAILURE);
    	  		}
    	  	}
    		}

        //Append seeds to allSeeds
        PList_appPList(allSeeds,seeds);
      }
    }

//  get info of growth curves
//  create an empty list (allGrowthVertices) for storing growth vertices of all growth curves
    pPList allGrowthVertices = PList_new();

    int ngc = PList_size(allSeeds);

    o.nGrowthCurves = ngc;
    o.arrays.gcflt = new double[ngc];
    o.arrays.gcgr  = new double[ngc];
    o.arrays.igcnv = new int[ngc];

    pPList growthVertices = PList_new();
    pPList growthEdges = PList_new();

//  FOR each seed in allSeeds
    for(int i = 0; i < PList_size(allSeeds); i++){
      seed = (pEdge)PList_item(allSeeds,i);

      PList_clear(growthVertices);
      PList_clear(growthEdges);

//    get growth vertices (growthVertices) and edges for seed
      if(!(BL_growthVerticesAndEdges((pEdge)seed, growthVertices, growthEdges) == 1)){
        printf("%s: unexpected BL_growthVerticesAndEdges return value\n",__func__);
        exit(EXIT_FAILURE);
      }

//    append growthVertices to allGrowthVertices
      PList_appPList(allGrowthVertices, growthVertices);

      o.arrays.igcnv[i] = PList_size(growthVertices);

      double l0 = E_length((pEdge)PList_item(growthEdges,0));
      o.arrays.gcflt[i] = l0;

      if( PList_size(growthEdges) > 1 )
        o.arrays.gcgr[i] = E_length((pEdge)PList_item(growthEdges,1))/l0;
      else
        o.arrays.gcgr[i] = 1.0;
    }

//  get info growth curves
    int nv = PList_size(allGrowthVertices);

    o.nLayeredMeshVertices = nv;
    o.arrays.igclv = new apf::MeshEntity*[nv];

    for(int i = 0; i < PList_size(allGrowthVertices); i++){
      vertex = (pVertex)PList_item(allGrowthVertices,i);

      apf::MeshEntity* me = reinterpret_cast<apf::MeshEntity*> (vertex);
      o.arrays.igclv[i] = me;
    }

    printf("%s: rank %d, ngc, nv: %d, %d\n", __func__, PCU_Comm_Self(), ngc, nv);

    PCU_Add_Ints(&ngc,sizeof(ngc));
    PCU_Add_Ints(&nv,sizeof(nv));

    if(PCU_Comm_Self() == 0)
      printf("%s: total ngc, nv: %d, %d\n", __func__, ngc, nv);

    PList_delete(gEdges);
    PList_delete(gVertices);
    PList_delete(gEntities);
    PList_delete(gFaces);
    PList_delete(seeds);
    PList_delete(allSeeds);
    PList_delete(blendSeeds);
    PList_delete(growthVertices);
    PList_delete(growthEdges);
    PList_delete(allGrowthVertices);

    //clean up utility
    Progress_delete(progress);
    Sim_logOff();
  }
  else {
    if(PCU_Comm_Self() == 0)
      printf("%s: warning! not implemented for MDS mesh\n",__func__);
  }
  return;
}
} //end namespace ph
