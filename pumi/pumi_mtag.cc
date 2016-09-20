/****************************************************************************** 

  (c) 2004-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "pumi.h"

pMeshTag pumi_mesh_createTag (pMesh mesh, const char* tagName, int tagType, int tagSize);
void pumi_mesh_deleteTag (pMesh mesh, pMeshTag tag, int forceDel);
pMeshTag pumi_mesh_findTag (pMesh mesh, const char* tagName);
bool pumi_mesh_hasTag (pMesh mesh, const pMeshTag tag);
void pumi_mesh_getTag (pMesh mesh, std::vector<pMeshTag>& tags);
bool pumi_mesh_isTagInUse (pMesh mesh, const pMeshTag tag); 
// sync tag data attached to the part boundary entities
void pumi_mesh_syncTag (pMesh mesh, pMeshTag tag, int ent_type);

void PUMI_MeshEnt_deleteTag (pMeshEnt meshEnt, pMeshTag tag);
bool PUMI_MeshEnt_hasTag (pMeshEnt meshEnt, pMeshTag tag);
void pumi_ment_getTag (pMeshEnt meshEnt, std::vector<pMeshTag>& tags);

//int pumi_ment_setPtrTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, void* data);
//int pumi_ment_getPtrTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, void** data);
void pumi_ment_setIntTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, const int data);
int pumi_ment_getIntTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, int* data);
void pumi_ment_setLongTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, const int data);
int pumi_ment_getLongTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, int* data);
void pumi_ment_setDblTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, const double data);
int pumi_ment_getDblTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, double* data);
//int pumi_ment_setEntTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, const pMeshEnt data);
//int pumi_ment_getEntTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, pMeshEnt* data);
//int pumi_ment_setSetTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, const pEntSet data);
//int pumi_ment_getSetTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, pEntSet* data);

//int pumi_ment_setPtrArrTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, void* const* data);
//int pumi_ment_getPtrArrTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, void** data);
void pumi_ment_setIntArrTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, const int* data, int data_size);
void pumi_ment_getIntArrTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, int** data, int* data_size);
void pumi_ment_setDblArrTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, const double* data, int data_size);
void pumi_ment_getDblArrTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, double** data, int* data_size);
//int pumi_ment_setEntArrTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, const pMeshEnt* data, int data_size);
//int pumi_ment_getEntArrTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, pMeshEnt** data, int* data_size);
//int pumi_ment_setSetArrTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, const pEntSet* data, int data_size);
//int pumi_ment_getSetArrTag (pMesh mesh, pMeshEnt ent, pMeshTag tag, pEntSet** data, int* data_size);
