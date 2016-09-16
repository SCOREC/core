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
#include "GenTag.h"
#include "GenIterator.h"
#include "mPartEntityContainer.h"

class gEntity;
class mPartEntityContainer;


class gModel : public TagHolder
{
private:  
  mPartEntityContainer allEntities;
  gmi_model* g; 
public:
  gModel(gmi_model* model);
  ~gModel();
  gmi_model* getGmi() {return g; }
  void add (int d, gEntity *ge) {allEntities.add(d, ge);}
  void del(int d, gEntity *ge) {allEntities.del(d, ge); } 
  typedef mPartEntityContainer::iter iterall;
  iterall beginall(int d) {return allEntities.begin(d);}
  iterall endall(int d) {return allEntities.end(d);}
  int size(int d) {return allEntities.size(d); }
};

typedef gModel* pGeom;
typedef gEntity* pGeomEnt;

typedef GenIterator<mPartEntityContainer::iter, gEntity>* gIter;

typedef apf::Mesh2* pMesh;
typedef apf::MeshEntity* pMeshEnt;
typedef apf::MeshEntity* pPartEnt;
typedef apf::MeshIterator* pMeshIter;
typedef apf::Copies Copies;
typedef apf::MeshTag* pMeshTag;
typedef apf::Parts Parts;
typedef apf::EntityVector EntityVector;
typedef apf::Parts Parts;
typedef apf::Up Up;
typedef apf::Downward Downward;
typedef apf::Migration Migration;

// singleton to save model/mesh
class pumi
{
public:
  pumi();
  ~pumi();
  static pumi* instance();

  pMesh mesh;
  pGeom model;
  pMeshTag ghosted_tag;
  pMeshTag ghost_tag;
  std::vector<pMeshEnt> ghost_vec[4];
  std::vector<pMeshEnt> ghosted_vec[4];
private:
  static pumi* _instance;
};

/** \brief Distribution plan object: send local elements to multiple destinations. */
// defined in pumi_distribution.cc
class Distribution
{
  public:
/** \brief must be constructed with a mesh
  \details use (new Distribution(mesh)) to make these objects */
    Distribution(pMesh m);
    ~Distribution();

/** \brief return true if the i'th element has been assigned destination(s) */
    bool has(pMeshEnt e);
/** \brief assign a destination part id to an element */
    void send(pMeshEnt e, int to);
/** \brief return the destination part id of an element */
    Parts& sending(pMeshEnt e);
    void print();
    int count();
    pMesh getMesh() {return m;}

    Parts* parts_vec;
    int element_count;
  private:
    pMesh m;
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
void pumi_printSys();
double pumi_getTime();
double pumi_getMem();
void pumi_printTimeMem(const char* msg, double time, double memory);

//************************************
//************************************
//      1- MESH FUNCTIONS
//************************************
//************************************

int pumi_tag_getType (const pTag tag);
void pumi_tag_getName (const pTag tag, const char** name);
int pumi_tag_getSize (const pTag tag);
int pumi_tag_getByte (const pTag tag);

//************************************
// Model management
//************************************

// create a model from a file
pGeom pumi_geom_load(const char* fileName, const char* model_type="mesh");

// iterator
/*
int pumi_giter_init (pGeom part, int type, gIter&);
int pumi_giter_getNext(gIter iter, pGeomEnt&);
void pumi_giter_delete(gIter iter);
void pumi_giter_reset(gIter iter);
bool pumi_giter_isEnd(gIter iter);
*/

pTag pumi_geom_createTag (pGeom g, const char* tagName, int tagType, int tagSize);
void pumi_geom_deleteTag (pGeom g, pTag tag, int forceDel=0);
pTag pumi_geom_findTag (pGeom g, const char* tagName);
bool pumi_geom_hasTag (pGeom g, const pTag tag);
bool pumi_geom_isTagInUse (pGeom g, const pTag tag);
void pumi_geom_getTag (pGeom g, std::vector<pTag>& tags);

bool pumi_gent_hasTag (pGeomEnt ent, pTag tag);
void pumi_gent_deleteTag (pGeomEnt ent, pTag tag);
void pumi_gent_getTag (pGeomEnt ent, std::vector<pTag>& tags);

//void pumi_gent_setStringTag(pGeomEnt ent, pTag tag, const char* s);
//void pumi_gent_getStringTag(pGeomEnt ent, pTag tag, const char*& s);

void pumi_gent_setPtrTag (pGeomEnt ent, pTag tag, void* data);
void pumi_gent_getPtrTag (pGeomEnt ent, pTag tag, void** data);
void pumi_gent_setIntTag (pGeomEnt ent, pTag tag, const int data);
void pumi_gent_getIntTag (pGeomEnt ent, pTag tag, int* data);
void pumi_gent_setLongTag (pGeomEnt ent, pTag tag, const long data);
void pumi_gent_getLongTag (pGeomEnt ent, pTag tag, long*);
void pumi_gent_setDblTag (pGeomEnt ent, pTag tag, const double data);
void pumi_gent_getDblTag (pGeomEnt ent, pTag tag, double*);
void pumi_gent_setEntTag (pGeomEnt ent, pTag tag, const pGeomEnt data);
void pumi_gent_getEntTag (pGeomEnt ent, pTag tag, pGeomEnt*);

void pumi_gent_setPtrArrTag (pGeomEnt ent, pTag tag, void* const* data);
void pumi_gent_getPtrArrTag (pGeomEnt ent, pTag tag, void** data);
void pumi_gent_setIntArrTag (pGeomEnt ent, pTag tag, const int* data);
void pumi_gent_getIntArrTag (pGeomEnt ent, pTag tag, int** data, int* data_size);
void pumi_gent_setLongArrTag (pGeomEnt ent, pTag tag, const long* data);
void pumi_gent_getLongArrTag (pGeomEnt ent, pTag tag, long** data, int* data_size);
void pumi_gent_setDblArrTag (pGeomEnt ent, pTag tag, const double* data);
void pumi_gent_getDblArrTag (pGeomEnt ent, pTag tag, double** data, int* data_size);
void pumi_gent_setEntArrTag (pGeomEnt ent, pTag tag, const pGeomEnt* data);
void pumi_gent_getEntArrTag (pGeomEnt ent, pTag tag, pGeomEnt** data, int* data_size);


//************************************
// Mesh management
//************************************

pGeom pumi_mesh_getGeom(pMesh m);

// load a serial mesh. 
pMesh pumi_mesh_loadSerial(pGeom g, const char* filename, const char* mesh_type="mds");

// load a mesh from a file. Do static partitioning if num_in_part==1
pMesh pumi_mesh_load(pGeom geom, const char* fileName, int num_in_part, const char* mesh_type="mds");

// load a serial mesh on master process then distribute as per the distribution object
void pumi_mesh_distribute(pMesh m, Distribution* plan);

// get mesh dimension
int pumi_mesh_getDim(pMesh m);

// get # mesh entities of type d on local process
int pumi_mesh_getNumEnt(pMesh m, int d);

// print mesh size info - global and local
void pumi_mesh_print(pMesh m, int p=0);

// write mesh into a file - mesh_type should be "mds" or "vtk"
void pumi_mesh_write (pMesh m, const char* fileName, const char* mesh_type="mds");

// delete mesh
void pumi_mesh_delete(pMesh m);

// verify mesh
void pumi_mesh_verify(pMesh m);

//************************************
// mesh tag management
//************************************

/*
pTag pumi_mesh_createTag (pMesh g, const char* tagName, int tagType, int tagSize);
void pumi_mesh_deleteTag (pMesh g, pTag tag, int forceDel=0);
pTag pumi_mesh_findTag (pMesh g, const char* tagName);
bool pumi_mesh_hasTag (pMesh g, const pTag tag);
bool pumi_mesh_isTagInUse (pMesh g, const pTag tag);
void pumi_mesh_getTag (pMesh g, std::vector<pTag>& tags);
// sync tag data attached to the part boundary entities
void pumi_mesh_syncTag (pMesh mesh, pMeshTag tag, int ent_type);

bool pumi_ment_hasTag (pMeshEnt ent, pTag tag);
void pumi_ment_deleteTag (pMeshEnt ent, pTag tag);
void pumi_ment_getTag (pMeshEnt ent, std::vector<pTag>& tags);

//void pumi_ment_setStringTag(pMeshEnt ent, pTag tag, const char* s);
//void pumi_ment_getStringTag(pMeshEnt ent, pTag tag, const char*& s);

void pumi_ment_setPtrTag (pMeshEnt ent, pTag tag, void* data);
void pumi_ment_getPtrTag (pMeshEnt ent, pTag tag, void** data);
void pumi_ment_setIntTag (pMeshEnt ent, pTag tag, const int data);
void pumi_ment_getIntTag (pMeshEnt ent, pTag tag, int* data);
void pumi_ment_setLongTag (pMeshEnt ent, pTag tag, const long data);
void pumi_ment_getLongTag (pMeshEnt ent, pTag tag, long*);
void pumi_ment_setDblTag (pMeshEnt ent, pTag tag, const double data);
void pumi_ment_getDblTag (pMeshEnt ent, pTag tag, double*);
void pumi_ment_setEntTag (pMeshEnt ent, pTag tag, const pMeshEnt data);
void pumi_ment_getEntTag (pMeshEnt ent, pTag tag, pMeshEnt*);

void pumi_ment_setPtrArrTag (pMeshEnt ent, pTag tag, void* const* data);
void pumi_ment_getPtrArrTag (pMeshEnt ent, pTag tag, void** data);
void pumi_ment_setIntArrTag (pMeshEnt ent, pTag tag, const int* data);
void pumi_ment_getIntArrTag (pMeshEnt ent, pTag tag, int** data, int* data_size);
void pumi_ment_setLongArrTag (pMeshEnt ent, pTag tag, const long* data);
void pumi_ment_getLongArrTag (pMeshEnt ent, pTag tag, long** data, int* data_size);
void pumi_ment_setDblArrTag (pMeshEnt ent, pTag tag, const double* data);
void pumi_ment_getDblArrTag (pMeshEnt ent, pTag tag, double** data, int* data_size);
void pumi_ment_setEntArrTag (pMeshEnt ent, pTag tag, const pMeshEnt* data);
void pumi_ment_getEntArrTag (pMeshEnt ent, pTag tag, pMeshEnt** data, int* data_size);
*/

//************************************
//  Ghosting
//************************************

/** \brief Ghosting plan object: local elements or part to destinations. */
// defined in pumi_ghost.cc
class Ghosting
{
  public:
/** \brief must be constructed with a mesh
  \details use (new apf::Migration(mesh)) to make these objects */
    Ghosting(pMesh, int d);
    ~Ghosting();

/** \brief return the number of elements with ghost destinations */
    int count();
    int count(pMeshEnt e, int d);
/** \brief return true if the i'th element has been assigned a destination */
    bool has(pMeshEnt e);
/** \brief assign a destination part id to an entity */
    void send(pMeshEnt e, int to);
/** \brief assign a destination part id of all ghost_dim entities */
    void send(int to);
    void print();
/** \brief return the destination parts of an entity */
    Parts& sending(pMeshEnt e, int d);

    pMesh getMesh() {return m;}
    int ghost_dim;

  private:
    pMesh m;
    pMeshTag parts_index_tag;
    std::vector<Parts*> parts_vec[4];
};

// migrate mesh per migration plan which contains a set of pairs [element and destination part]
void pumi_mesh_migrate(pMesh m, Migration* plan);

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
void pumi_ghost_createLayer (pMesh m, int brgType, int ghostType, int numLayer, int includeCopy);

/* 
  Ghosting: ghosting plan object for local elements or part to destinations.
*/
void pumi_ghost_create(pMesh m, Ghosting* plan);

// unavailable
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
// unavailable
void pumi_ghost_getInfo (pMesh m, std::vector<int>& ghostinfo);

//************************************
//  Mesh Entity
//************************************
// get geometric entity's dimension
int pumi_gent_getDim(pGeomEnt ge);

// get geometric entity's global id
int pumi_gent_getID(pGeomEnt ge);

void pumi_gent_getRevClas (pGeomEnt g, std::vector<pMeshEnt>& ents);

// get mesh entity's dimension
int pumi_ment_getDim(pMeshEnt e);

// get mesh entity's local id
int pumi_ment_getLocalID(pMeshEnt e);

// get mesh entity's global id - vertex only
// global id is maintained if mesh is ghosted
// global id is NOT maintained if mesh is adapted or re-partitioned
int pumi_ment_getGlobalID(pMeshEnt e);

// get # adjacent entities
int pumi_ment_getNumAdj(pMeshEnt e, int tgtType);

// get adjacent entities
void pumi_ment_getAdj(pMeshEnt e, int tgtType, std::vector<pMeshEnt>& vecAdjEnt);

// get 2nd-order adjacent entities
void pumi_ment_get2ndAdj (pMeshEnt e, int brgType, int tgtType, std::vector<pMeshEnt>& vecAdjEnt);

// return entity's geometric classification
pGeomEnt pumi_ment_getGeomClas(pMeshEnt e);

// unavailable
pPartEnt pumi_ment_getPtnClas(pMeshEnt e);

// return owning part id. if ghosted mesh, vertex or element only
int pumi_ment_getOwnPID(pMeshEnt e); 

// return owner entity copy. if ghoted mesh, vertex or element only
pMeshEnt pumi_ment_getWwnEnt(pMeshEnt e); 

// return true if the entity is an owner copy
bool pumi_ment_isOwned(pMeshEnt e);

// return true if entity is on part boundary, ghosted, or ghost
//  - this will fixed to consider only part boundary entities later
bool pumi_ment_isOnBdry (pMeshEnt e); 

// return # remote and ghost copies
//  - this will fixed to consider only part boundary entities later
int pumi_ment_getNumRmt (pMeshEnt e); 

// return remote and ghost copies
//  - this will fixed to consider only part boundary entities later
void pumi_ment_getAllRmt(pMeshEnt e, Copies& remotes); 

// return remote or ghost copy on a destination part
//  - this will fixed to consider only part boundary entities later
pMeshEnt pumi_ment_getRmt(pMeshEnt& meshEnt, int destPart); 

// unavailable
void pumi_ment_setRmt(pMeshEnt e, int partID, pMeshEnt rmtEnt);

// unavailable
void pumi_ment_deleteRmt (pMeshEnt e, int partID);

// unavailable
void pumi_ment_cleanRmt (pMeshEnt e);

// unavailable
void pumi_ment_setPtnTopology (pMeshEnt e);

// return part ids where the entity is duplicated - part boundary or ghost
void pumi_ment_getResidence(pMeshEnt e, std::vector<int>& resPartId);

// return part ids where the entity and its downward adjacency are duplicated - part boundary or ghost
void pumi_ment_getClosureResidence(pMeshEnt ent, std::vector<int>& resPartId);

// return true if the entity is a ghost copy
bool pumi_ment_isGhost(pMeshEnt e);

// return true if the entity is ghosted
// unavailable
bool pumi_ment_isGhosted (pMeshEnt e);

// unavailable
// return #ghost copies
int pumi_ment_getNumGhost (pMeshEnt e);

// unavailable
// return ghost copies
void pumi_ment_getAllGhost (pMeshEnt e, Copies&);

// unavailable
// return ghost copy on a destination part
pMeshEnt pumi_ment_getGhost(pMeshEnt& e, int partID);
#endif
