/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef PUMI_H
#define PUMI_H
#include <gmi.h>
#include <apfMesh2.h>

typedef gmi_model* pGeom;
typedef gmi_ent* pGeomEnt;
typedef apf::Mesh2* pMesh;
typedef apf::MeshEntity* pMeshEnt;
typedef apf::MeshEntity* pPartEnt;
typedef apf::MeshIterator* pMeshIter;
typedef apf::Copies pCopies;
typedef apf::MeshTag* pTag;

// singleton to save model/mesh
class pumi
{
public:
  pumi();
  ~pumi();
  static pumi* instance();

  pMesh mesh;
  pGeom model;
  std::map<pMeshEnt, pCopies> ghost_map[4];
  std::map<pMeshEnt, std::set<int> > bps_map[4];
  pMesh org_mesh;
  std::vector<bool>* org_node_flag;
  int num_own_vtx;
private:
  static pumi* _instance;
};

//************************************
//************************************
//      0- SYSTEM-LEVEL FUNCTIONS
//************************************
//************************************
void pumi_start();
void pumi_finalize(bool do_mpi_finalize=true);

int pumi_size();
int pumi_rank();

void pumi_sync(void);
void pumi_printsys();
double pumi_gettime();
double pumi_getmem();
void pumi_printtimemem(const char* msg, double time, double memory);

//************************************
//************************************
//      1- MESH FUNCTIONS
//************************************
//************************************

//************************************
// Model/Mesh management
//************************************
// create a model from a file
pGeom pumi_geom_create(const char* fileName, const char* model_type="mesh");

// create a mesh from a file
pMesh pumi_mesh_create(pGeom geom, const char* fileName, int option, int num_proc_grp=1, const char* mesh_type="mds");

// get mesh dimension
int pumi_mesh_getdim(pMesh m);

// print mesh size info - global and local
void pumi_mesh_print(pMesh m);

// write mesh into a file
void pumi_mesh_write (pMesh m, const char* fileName, const char* mesh_type="mds");

// delete mesh
void pumi_mesh_delete(pMesh m);

// get # local entities of type d
int pumi_mesh_getnument(pMesh m, int d);

//************************************
//  Ghosting
//************************************
/* 
input:
  - brgType - desired bridge entity type
  - ghostType - desired ghost entity type 
  - numLayer - the number of ghost layers
  - includeCopy - integer flag indicating whether to include non-owned bridge entity (1: yes, 0: no)
	If includeCopy is 0 and part boundary entity of brgType is not owned by a self part 
	(shortly, non-owned bridge type entity), ghost type dimensional entities adjacent to the non-owned
	bridge type entity is not ghost copied. If includeCopy is 1, all ghost type dimensional entities 
	adjacent to the bridge type entities on part boundaries are ghost copied.

The error is returned in the following cases:
  - bridge type is greater than or equal to ghost type
  - bridge type is greater than or equal to mesh dimension
  - ghost type is mesh vertex
  - ghost type is grester than mesh dimension
*/
pMesh pumi_ghost_create (int brgType, int ghostType, int numLayer, int includeCopy);

void pumi_ghost_delete (pMesh m);

/* 
input: 	
  - a mesh instance
   
output: 
  - return the historical ghosting information in order which consists of four integers, ghost type, 
    bridge type, the number of layers, and includeCopy flag (0 or 1)

example:
  If pumi_ghost_create was called twice in the following order (mesh, 0, 2, 1, 1) and (mesh, 1, 3, 0), 
  the vector "ghostinfo" contains [0, 2, 1, 1, 1, 3, 3, 0]
*/
void pumi_ghost_info (pMesh m, std::vector<int>& ghostinfo);

//************************************
//  Mesh Entity
//************************************
int pumi_gent_getdim(pGeomEnt ge);
int pumi_gent_getid(pGeomEnt ge);

int pumi_ment_getdim(pMeshEnt e);
int pumi_ment_getid(pMeshEnt e);

int pumi_ment_getnumadj(pMeshEnt e, int tgtType);
void pumi_ment_getadj(pMeshEnt e, int tgtType, std::vector<pMeshEnt>& vecAdjEnt);
void pumi_ment_get2ndadj (pMeshEnt e, int brgType, int tgtType, std::vector<pMeshEnt>& vecAdjEnt);

pGeomEnt pumi_ment_getgeomclass(pMeshEnt e);
pPartEnt pumi_ment_getptnclass(pMeshEnt e);

int pumi_ment_getglobalid(pMeshEnt e); // vertex only

// owner part information
int pumi_ment_getownpid(pMeshEnt e);
pMeshEnt pumi_ment_getownent(pMeshEnt e);
bool pumi_ment_isowned(pMeshEnt e);

// remote copy information
bool pumi_ment_isonbdry (pMeshEnt e);
int pumi_ment_getnumrmt (pMeshEnt e);
void pumi_ment_getallrmt(pMeshEnt e, pCopies& remotes);
pMeshEnt pumi_ment_getrmt(pMeshEnt& meshEnt, int destPart);
void pumi_ment_setrmt(pMeshEnt e, int partID, pMeshEnt rmtEnt);
void pumi_ment_deletermt (pMeshEnt e, int partID);
void pumi_ment_cleanrmt (pMeshEnt e);

void pumi_ment_setptntopology (pMeshEnt e);
void pumi_ment_getresidence(pMeshEnt e, std::vector<int>& resPartId);
void pumi_ment_getclosureresidence(pMeshEnt ent, std::vector<int>& resPartId);

// ghosting information
bool pumi_ment_isghost(pMeshEnt e);
bool pumi_ment_isghosted (pMeshEnt e);
int pumi_ment_getnumghost (pMeshEnt e);
void pumi_ment_getallghost (pMeshEnt e, pCopies&);
pMeshEnt pumi_ment_getghost(pMeshEnt& e, int partID);
#endif
