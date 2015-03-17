/*
 * Copyright (C) 2014 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfZoltanCallbacks.h"
#include "apfZoltanMesh.h"
#include "apfZoltan.h"
#include "apfShape.h"
#include <PCU.h>
#include <metis.h>

namespace apf {

static int setZoltanLbMethod(struct Zoltan_Struct* ztn, ZoltanMesh* zb)
{
  // setting LB_METHOD
  std::string lbMethod = "GRAPH";
  switch (zb->method) {
    case RCB:
      lbMethod = "RCB"; break;
    case RIB:
      lbMethod = "RIB"; break;
    case HYPERGRAPH:
      lbMethod = "HYPERGRAPH"; break;
    case PARMETIS: //fall into GRAPH settings
    case GRAPH:
      lbMethod = "GRAPH";
      Zoltan_Set_Param(ztn, "GRAPH_PACKAGE", "PARMETIS"); // instead of PHG
      break;
    default:
      std::cout << "ERROR " << __func__ << " Invalid LB_METHOD "
        << zb->method << "\n";
      return 1;
  }
  Zoltan_Set_Param(ztn, "LB_METHOD", lbMethod.c_str());
  return 0;
}

static int setZoltanLbApproach(struct Zoltan_Struct* ztn, ZoltanMesh* zb)
{
  // setting LB_Approach
  std::string ptnAp = "REPARTITION";
  std::string pMethod = "PartKway";
  switch (zb->approach) {
    case PARTITION:
      ptnAp = "PARTITION"; break;
    case REPARTITION:
      ptnAp = "REPARTITION"; break;
    case REFINE:
      ptnAp = "REFINE"; break;
    case PART_KWAY:
      pMethod = "PartKway"; break;
    case PART_GEOM:
      pMethod = "PartGeom"; break;
    case PART_GEOM_KWAY:
      pMethod = "PartGeomKway"; break;
    case ADAPT_REPART:
      pMethod = "AdaptiveRepart"; break;
    case REFINE_KWAY:
      pMethod = "RefineKway"; break;
    default:
      std::cout << "ERROR " << __func__ << " Invalid LB_Approach "
        << zb->approach << "\n";
      return 1;
  }
  Zoltan_Set_Param(ztn, "LB_APPROACH", ptnAp.c_str());
  if ( (3 == zb->method) || (4 == zb->method) )
    Zoltan_Set_Param(ztn, "PARMETIS_METHOD", pMethod.c_str());
  return 0;
}

static long get(MeshEntity* e, ZoltanMesh *zz)
{
  if (!zz->isLocal)
    return getNumber(zz->global,Node(e,0));
  else
    return getNumber(zz->local, e, 0, 0);
}

static int getPartId(Mesh* m, MeshEntity* s)
{
  if (m->isShared(s))
    return getOtherCopy(m, s).peer;
  Matches matches;
  m->getMatches(s, matches);
  assert(matches.getSize() == 1);
  return matches[0].peer;
}

//ZOLTAN_NUM_OBJ_FN_TYPE
int zoltanCountNodes(void* data, int* ierr)
{
  ZoltanMesh* zb = static_cast<ZoltanMesh*>(data);
  *ierr=ZOLTAN_OK;
  return zb->elements.getSize();
}

//ZOLTAN_OBJ_LIST_FN_TYPE
void zoltanGetNodes(void* data, int ngid, int nlid,
    ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int nweights,
    float* weights, int* ierr)
{
  ZoltanMesh* zb = static_cast<ZoltanMesh*>(data);
  *ierr=ZOLTAN_OK;
  for (size_t ind=0;ind<zb->elements.getSize();ind++) {
    lids[ind*nlid]=ind;
    MeshEntity *ent = zb->elements[ind];
    gids[ind*ngid]= get(ent, zb);
    double* w = (double*)calloc(nweights,sizeof(double));
    zb->mesh->getDoubleTag(ent,zb->weights,w);
    for (int i=0;i<nweights;i++)
      weights[ind*nweights+i]=w[i];
    free(w);
  }

  //if weights arent being used
  if (nweights==0)
    weights=NULL;
}

//ZOLTAN_NUM_EDGES_FN_TYPE
int zoltanCountEdges(void* data, int, int,
    ZOLTAN_ID_PTR, ZOLTAN_ID_PTR lid, int* ierr)
{
  ZoltanMesh* zb = static_cast<ZoltanMesh*>(data);
  int edges = 0;
  MeshEntity* element=zb->elements[*lid];
  Mesh* mesh = zb->mesh;
  int dim = getDimension(mesh, element);
  Downward faces;
  int nfaces;
  nfaces = mesh->getDownward(element,dim-1,faces);
  for (int i=0;i<nfaces;i++) {
    MeshEntity* face = faces[i];
    Up elements;
    mesh->getUp(face,elements);
    if (elements.n==2)
      edges++;
    else if (elements.n==1&&!zb->isLocal&&mesh->hasTag(face,zb->opposite))
      edges++;
  }
  *ierr = ZOLTAN_OK;
  return edges;
}

//ZOLTAN_EDGE_LIST_FN_TYPE
void zoltanGetEdges(void* data, int ngid, int,
    ZOLTAN_ID_PTR, ZOLTAN_ID_PTR lid,
    ZOLTAN_ID_PTR gids, int* pids,
    int, float*, int* ierr)
{
  ZoltanMesh* zb = static_cast<ZoltanMesh*>(data);
  MeshEntity* element=zb->elements[*lid];
  Mesh* mesh = zb->mesh;
  int dim = getDimension(mesh, element);
  Downward faces;
  int nfaces;
  nfaces = mesh->getDownward(element,dim-1,faces);
  int ind=0;
  for (int i=0;i<nfaces;i++) {
    MeshEntity* face = faces[i];
    Up elements;
    mesh->getUp(face,elements);
    if (elements.n==2) {
      MeshEntity* ent = elements.e[0];
      if (element==ent)
        ent=elements.e[1];
      gids[ngid*ind] = get(ent,zb);
      pids[ind] = zb->isLocal ? 0 : mesh->getId();
      ind++;
    }
    else if (elements.n==1&&!zb->isLocal&&mesh->hasTag(face,zb->opposite)) {
      long value;
      mesh->getLongTag(face,zb->opposite,&value);
      gids[ngid*ind]=value;
      pids[ind]=getPartId(mesh,face);
      ind++;
    }
  }
  *ierr = ZOLTAN_OK;
}

// ZOLTAN_GEOM_FN
void getCentroid(void *data, int, int,
    ZOLTAN_ID_PTR, ZOLTAN_ID_PTR lid, double *coords,
    int *ierr)
{
  ZoltanMesh* zb = static_cast<ZoltanMesh*>(data);
  getLinearCentroid(zb->mesh,zb->elements[*lid]).toArray(coords);
  *ierr=ZOLTAN_OK;
}

// ZOLTAN_NUM_GEOM_FN_TYPE
int getGeomDim(void *, int *)
{
  return 3; //always 3D!!
}

// ZOLTAN_HG_CS_FN
void getHg(void* data, int ngid, int, int totAdjVtx, int,
    ZOLTAN_ID_PTR elmIds, int* adjVtxIdx, ZOLTAN_ID_PTR adjVtx,
    int *ierr)
{
  ZoltanMesh* zb = static_cast<ZoltanMesh*>(data);
  Numbering* ln = NULL;
  GlobalNumbering* gn = NULL;
  if( zb->isLocal ) {
    ln = numberOverlapNodes(zb->mesh, "zoltan_vtx");
  } else {
    gn = makeGlobal(numberOwnedNodes(zb->mesh, "zoltan_vtx"));
    synchronize(gn);
  }
  Downward verts;
  int prevCnt = adjVtxIdx[0] = 0;
  for (size_t ind=0;ind<zb->elements.getSize();ind++) {
    MeshEntity* elm = zb->elements[ind];
    elmIds[ind*ngid]= get(elm, zb);
    adjVtxIdx[ind] = prevCnt;
    int nv = prevCnt = zb->mesh->getDownward(elm, 0, verts);
    prevCnt += adjVtxIdx[ind];
    for(int i=0; i<nv; i++) {
      //add get(...) support for vertices???
      const int avIdx = (adjVtxIdx[ind]+i)*ngid;
      assert(totAdjVtx > avIdx);
      long gid = zb->isLocal ?
        getNumber(ln,verts[i],0,0) : getNumber(gn,Node(verts[i],0));
      adjVtx[avIdx] = gid;
    }
  }
  if( zb->isLocal )
    destroyNumbering(ln);
  else
    destroyGlobalNumbering(gn);
  *ierr=ZOLTAN_OK;
}

// ZOLTAN_HG_SIZE_CS_FN
void getHgSize(void* data,
    int* numElms, //out
    int* numAdjVtx, //out - sum(getNumVtx(elm[i]))
    int* format, //out - hardcoded to compressed vtx
    int* ierr) {
  ZoltanMesh* zb = static_cast<ZoltanMesh*>(data);

  *format = ZOLTAN_COMPRESSED_VERTEX;
  *numElms = zb->elements.getSize();

  *numAdjVtx = 0;
  for (size_t i=0;i<zb->elements.getSize();i++) {
    MeshEntity* elm = zb->elements[i];
    const int type = zb->mesh->getType(elm);
    const int nsides = apf::Mesh::adjacentCount[type][0];
    *numAdjVtx += nsides;
  }
  *ierr=ZOLTAN_OK;
}

ZoltanData::ZoltanData(ZoltanMesh* zb_) : zb(zb_)
{
  // allocate zoltan data structures
  float ver;
  int ret = Zoltan_Initialize(0,0,&ver);
  if (ZOLTAN_OK != ret) {
    fprintf(stderr, "ERROR: Zoltan initialization failed\n");
    exit(1);
  }
  MPI_Comm comm;
  if (zb->isLocal)
    comm = MPI_COMM_SELF;
  else
    comm = PCU_Get_Comm();
  ztn =Zoltan_Create(comm);
  import_gids = NULL;
  import_lids = NULL;
  import_procs = NULL;
  export_gids = NULL;
  export_lids = NULL;
  export_procs = NULL;
  num_imported = 0;
  num_exported = 0;
  import_to_part = NULL;
  export_to_part = NULL;
  dbgLvl = zb->debug ? (1) : (0);
  changes=0;
  lidSz=1;
  gidSz=1;
}

ZoltanData::~ZoltanData()
{
  // delete zoltan data structures
  Zoltan_LB_Free_Part(&import_gids, &import_lids, &import_procs, &import_to_part);
  Zoltan_LB_Free_Part(&export_gids, &export_lids, &export_procs, &export_to_part);
  Zoltan_Destroy(&ztn);
}

void ZoltanData::run()
{
/* ensure Metis indices are the same size as Zoltan indices */
  assert(IDXTYPEWIDTH == sizeof(ZOLTAN_ID_TYPE)*8);
  setup();
  ptn();
}

void ZoltanData::setup()
{
  //fill out zoltan parameters
  char paramStr[128];

  //sizes
  sprintf(paramStr, "%u", 1);
  Zoltan_Set_Param(ztn, "num_gid_entries", paramStr);
  sprintf(paramStr, "%u", 1);
  Zoltan_Set_Param(ztn, "num_lid_entries", paramStr);

  //weights
  sprintf(paramStr, "%d", zb->mesh->getTagSize(zb->weights));
  Zoltan_Set_Param(ztn, "obj_weight_dim", paramStr);
  Zoltan_Set_Param(ztn, "edge_weight_dim", "0");

  //Debug
  sprintf(paramStr, "%d", dbgLvl);
  if ( zb->isLocal && 0 != PCU_Comm_Self() )
    sprintf(paramStr, "%d", 0);  //if local silence all but rank 0
  Zoltan_Set_Param(ztn, "debug_level", paramStr);
  Zoltan_Set_Param(ztn, "PARMETIS_OUTPUT_LEVEL", paramStr);
  Zoltan_Set_Param(ztn, "CHECK_GRAPH", "0");
  Zoltan_Set_Param(ztn, "CHECK_HYPERGRAPH", "0");

  //tolerance
  sprintf(paramStr, "%f", zb->tolerance);
  Zoltan_Set_Param(ztn, "imbalance_tol", paramStr);

  Zoltan_Set_Param(ztn, "RETURN_LISTS", "EXPORT");

  setZoltanLbMethod(ztn,zb);
  setZoltanLbApproach(ztn,zb);

  /* Reset some load-balancing parameters. */
  if ( zb->isLocal ) {
    sprintf(paramStr, "%d", zb->multiple);
  } else {
    sprintf(paramStr, "%d", zb->multiple*PCU_Proc_Peers());
  }
  Zoltan_Set_Param(ztn, "NUM_GLOBAL_PARTS", paramStr);
  sprintf(paramStr, "%d", zb->multiple);
  Zoltan_Set_Param(ztn, "NUM_LOCAL_PARTS", paramStr);

  Zoltan_Set_Param(ztn, "GRAPH_BUILD_TYPE", "FAST_NO_DUP");

  //set zoltan call backs
  Zoltan_Set_Fn(ztn, ZOLTAN_NUM_OBJ_FN_TYPE, (void (*)())zoltanCountNodes, (void*) (zb));
  Zoltan_Set_Fn(ztn, ZOLTAN_OBJ_LIST_FN_TYPE, (void (*)())zoltanGetNodes, (void *) (zb));
  Zoltan_Set_Fn(ztn, ZOLTAN_NUM_EDGES_FN_TYPE, (void (*)())zoltanCountEdges, (void*) (zb));
  Zoltan_Set_Fn(ztn, ZOLTAN_EDGE_LIST_FN_TYPE, (void (*)())zoltanGetEdges, (void*) (zb));
  Zoltan_Set_Fn(ztn, ZOLTAN_NUM_GEOM_FN_TYPE, (void (*)()) getGeomDim, (void*) (zb));
  Zoltan_Set_Fn(ztn, ZOLTAN_GEOM_FN_TYPE, (void (*)()) getCentroid, (void*) (zb));
  Zoltan_Set_Fn(ztn, ZOLTAN_HG_SIZE_CS_FN_TYPE, (void (*)()) getHgSize, (void*) (zb));
  Zoltan_Set_Fn(ztn, ZOLTAN_HG_CS_FN_TYPE, (void (*)()) getHg, (void*) (zb));
}

void ZoltanData::ptn()
{
  //run partitioner
  int ret = Zoltan_LB_Partition(ztn, &changes, &gidSz, &lidSz,
      &num_imported, &import_gids, &import_lids, &import_procs,
      &import_to_part, &num_exported, &export_gids,
      &export_lids, &export_procs, &export_to_part);
  if( ZOLTAN_OK != ret ) {
    if( 0 == PCU_Comm_Self() )
      fprintf(stderr, "ERROR Zoltan partitioning failed\n");
    exit(EXIT_FAILURE);
  }
  //writeZoltanDbgFiles(ztn, "kddParted");
}

void ZoltanData::getExport(int ind, int *localId, int *export_part)
{
  assert(export_lids[ind]<zb->elements.getSize());
  *localId = export_lids[ind];
  *export_part = export_to_part[ind];
}

}
