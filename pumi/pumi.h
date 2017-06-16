/****************************************************************************** 

  (c) 2004-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef PUMI_H
#define PUMI_H
#include <gmi.h>
#include <apfMesh2.h>
#include "GenTag.h"
#include "GenIterator.h"
#include "mPartEntityContainer.h"
#include "apf.h"

enum PUMI_EntTopology {
  PUMI_VERTEX, // 0 
  PUMI_EDGE,   // 1 
  PUMI_TRIANGLE, // 2 
  PUMI_QUAD, // 3
  PUMI_TET,  // 4
  PUMI_HEX,  // 5 
  PUMI_PRISM, // 6
  PUMI_PYRAMID, // 7
  ENT_TOPOLOGIES
};

enum PUMI_FieldType {
  PUMI_SCALAR, //  a single scalar value
  PUMI_VECTOR, // a 3D vector
  PUMI_MATRIX, // a 3x3 matrix 
  PUMI_PACKED, // a user-defined set of components
  FIELD_TYPES
};
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
  gEntity* getGeomEnt(int d, gmi_ent* ge);
  void add (int d, gEntity *ge) {allEntities.add(d, ge);}
  void del(int d, gEntity *ge) {allEntities.del(d, ge); } 
  typedef mPartEntityContainer::iter iterall;
  iterall begin(int d) {return allEntities.begin(d);}
  iterall end(int d) {return allEntities.end(d);}
  int size(int d) {return allEntities.size(d); }
};

typedef gModel* pGeom;
typedef gEntity* pGeomEnt;

typedef gModel::iterall pGeomIter;
typedef GenIterator<mPartEntityContainer::iter, gEntity>* gIter;

typedef apf::Mesh2* pMesh;
typedef apf::MeshEntity* pMeshEnt;
typedef apf::EntityVector EntityVector;
typedef apf::MeshIterator* pMeshIter;
typedef apf::Copies Copies;
typedef apf::Copies::iterator pCopyIter;
typedef apf::MeshTag* pMeshTag;
typedef apf::Parts Parts;
typedef apf::Migration Migration;
typedef apf::Field* pField;
typedef apf::FieldShape* pShape;
typedef apf::Numbering* pNumbering;
typedef apf::GlobalNumbering* pGlobalNumbering;
typedef apf::Vector3 Vector3; // 3d vector

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
//      1- MODEL FUNCTIONS
//************************************
//************************************

// Geometric Model
// create a model from a file
pGeom pumi_geom_load (const char* fileName, const char* model_type="mesh", 
                      void (*fp)(const char*)=NULL);
void pumi_geom_freeze(pGeom g); // shall be called after modifying model entities
int pumi_geom_getNumEnt(pGeom g, int d);
pGeomEnt pumi_geom_findEnt(pGeom g, int d, int id);

// Geometric Entity
// get geometric entity's dimension
int pumi_gent_getDim(pGeomEnt ge);
// get geometric entity's local id (note local==global for geometric model)
int pumi_gent_getID(pGeomEnt ge);
void pumi_gent_getRevClas (pGeomEnt g, std::vector<pMeshEnt>& ents);
int pumi_gent_getNumAdj (pGeomEnt g, int target_dim);
void pumi_gent_getAdj (pGeomEnt g, int target_dim, std::vector<pGeomEnt>& ents);

// Tag management
pTag pumi_geom_createTag (pGeom g, const char* tagName, int tagType, int tagSize);
void pumi_geom_deleteTag (pGeom g, pTag tag, bool force_delete=false);
pTag pumi_geom_findTag (pGeom g, const char* tagName);
bool pumi_geom_hasTag (pGeom g, const pTag tag);
void pumi_geom_getTag (pGeom g, std::vector<pTag>& tags);

int pumi_tag_getType (const pTag tag);
void pumi_tag_getName (const pTag tag, const char** name);
int pumi_tag_getSize (const pTag tag);
int pumi_tag_getByte (const pTag tag);

void pumi_gent_deleteTag (pGeomEnt ent, pTag tag);
bool pumi_gent_hasTag (pGeomEnt ent, pTag tag);
void pumi_gent_getTag (pGeomEnt ent, std::vector<pTag>& tags);

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

// create an empty mesh
pMesh pumi_mesh_create(pGeom g, int mesh_dim, bool periodic=false);
void pumi_mesh_freeze(pMesh m);
pMeshEnt pumi_mesh_createVtx(pMesh m, pGeomEnt ge, double* xyz);
//ent_topology: VERTEX (0), EDGE (1), TRIANGLE (2), QUAD (3), TET (4), HEX (5), PRISM (6), PYRAMID (7)
pMeshEnt pumi_mesh_createEnt(pMesh m, pGeomEnt ge, int target_topology, pMeshEnt* down);
pMeshEnt pumi_mesh_createElem(pMesh m, pGeomEnt ge, int target_topology, pMeshEnt* vertices); 

// load a serial mesh. 
pMesh pumi_mesh_loadSerial(pGeom g, const char* file_name, const char* mesh_type="mds");

// load a mesh from a file. Do static partitioning if num_in_part==1
pMesh pumi_mesh_load(pGeom geom, const char* fileName, int num_in_part, const char* mesh_type="mds");

// delete mesh
void pumi_mesh_delete(pMesh m);
// write mesh into a file - mesh_type should be "mds" or "vtk"
void pumi_mesh_write (pMesh m, const char* fileName, const char* mesh_type="mds");
pGeom pumi_mesh_getGeom(pMesh m);

// get mesh dimension
int pumi_mesh_getDim(pMesh m);

// get # mesh entities of type d on local process
int pumi_mesh_getNumEnt(pMesh m, int d);
pMeshEnt pumi_mesh_findEnt(pMesh m, int d, int id);

//************************************
// Tag management over mesh
//************************************
pMeshTag pumi_mesh_createIntTag(pMesh m, const char* name, int size);
pMeshTag pumi_mesh_createLongTag(pMesh m, const char* name, int size);
pMeshTag pumi_mesh_createDblTag(pMesh m, const char* name, int size);

void pumi_mesh_deleteTag(pMesh m, pMeshTag tag, bool force_delete=false);
pMeshTag pumi_mesh_findTag(pMesh m, const char* name);
bool pumi_mesh_hasTag (pMesh m, const pMeshTag tag);
void pumi_mesh_getTag(pMesh m, std::vector<pMeshTag> tags);

//************************************
//  Migration
//************************************

// migrate mesh per migration plan which contains a set of pairs [element and destination part]
void pumi_mesh_migrate(pMesh m, Migration* plan);

//************************************
//  Distribution
//************************************

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

// load a serial mesh on master process then distribute as per the distribution object
void pumi_mesh_distribute(pMesh m, Distribution* plan);

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

/* 
input:
  - brgType - desired bridge entity dimension
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

// Ghosting: ghosting plan object for local elements or part to destinations. 
void pumi_ghost_create(pMesh m, Ghosting* plan);

void pumi_ghost_delete (pMesh m);

//************************************
// MISCELLANEOUS
//************************************
// create/delete global ID using tag "global_id"
void pumi_mesh_createGlobalID(pMesh m);
void pumi_mesh_deleteGlobalID(pMesh m);

// verify mesh
void pumi_mesh_verify(pMesh m, bool abort_on_error=true);

// print mesh size info - global and local
void pumi_mesh_print(pMesh m, int print_ent=0);

//************************************
//  Mesh Entity
//************************************
// get mesh entity's dimension
int pumi_ment_getDim(pMeshEnt e);
int pumi_ment_getTopo(pMeshEnt e);
// get mesh entity's local id
int pumi_ment_getID(pMeshEnt e);

// get mesh entity's global id - vertex only
// global id is maintained if mesh is re-partitioned or ghosted
// global id is NOT maintained if mesh is adapted
int pumi_ment_getGlobalID(pMeshEnt e);

// get # adjacent entities
int pumi_ment_getNumAdj(pMeshEnt e, int tgtType);

// get adjacent entities
void pumi_ment_getAdj(pMeshEnt e, int tgtType, std::vector<pMeshEnt>& vecAdjEnt);

// get 2nd-order adjacent entities
void pumi_ment_get2ndAdj (pMeshEnt e, int brgType, int tgtType, std::vector<pMeshEnt>& vecAdjEnt);

// return entity's geometric classification
pGeomEnt pumi_ment_getGeomClas(pMeshEnt e);

// for mesh edge and vertex, return the other vertex
pMeshEnt pumi_medge_getOtherVtx(pMeshEnt edge, pMeshEnt vtx);

// tag management over mesh entity
void pumi_ment_deleteTag (pMeshEnt e, pMeshTag tag);
bool pumi_ment_hasTag (pMeshEnt e, pMeshTag tag);

void pumi_ment_setIntTag(pMeshEnt e, pMeshTag tag, int const* data);
void pumi_ment_getIntTag(pMeshEnt e, pMeshTag tag, int* data);
void pumi_ment_setLongTag(pMeshEnt e, pMeshTag tag, long const* data);
void pumi_ment_getLongTag(pMeshEnt e, pMeshTag tag, long* data);
void pumi_ment_setDblTag(pMeshEnt e, pMeshTag tag, double const* data);
void pumi_ment_getDblTag(pMeshEnt e, pMeshTag tag, double* data);

// return owning part id. if ghosted mesh, vertex or element only
int pumi_ment_getOwnPID(pMeshEnt e); 

// return owner entity copy. if ghoted mesh, vertex or element only
pMeshEnt pumi_ment_getOwnEnt(pMeshEnt e); 

// return true if the entity is an owner copy
bool pumi_ment_isOwned(pMeshEnt e);

// return true if the entity exists on the part
bool pumi_ment_isOn(pMeshEnt e, int partID);

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

// return part ids where the entity is duplicated - part boundary or ghost
void pumi_ment_getResidence(pMeshEnt e, Parts& residence);

// return part ids where the entity and its downward adjacency are duplicated - part boundary or ghost
void pumi_ment_getClosureResidence(pMeshEnt ent, Parts& residence);

// return true if the entity is a ghost copy
bool pumi_ment_isGhost(pMeshEnt e);

// return true if the entity is ghosted
bool pumi_ment_isGhosted (pMeshEnt e);

// return #ghost copies
int pumi_ment_getNumGhost (pMeshEnt e);

// return ghost copies
void pumi_ment_getAllGhost (pMeshEnt e, Copies&);

// return ghost copy on a destination part
pMeshEnt pumi_ment_getGhost(pMeshEnt& e, int partID);

//************************************
// Mesh entity numbering
//************************************

void pumi_ment_setGlobalNumber(pMeshEnt e, pGlobalNumbering gn, int node, int component, long number);
long pumi_ment_getGlobalNumber(pMeshEnt e, pGlobalNumbering gn, int node=0, int component=0);
void pumi_ment_setNumber(pMeshEnt e, pNumbering n, int node, int component, int number);
int pumi_ment_getNumber(pMeshEnt e, pNumbering n, int node=0, int component=0);
bool pumi_ment_isNumbered(pMeshEnt e, pNumbering n);

//************************************
// Field shape and nodes
//************************************
pShape pumi_mesh_getShape (pMesh m);
void pumi_mesh_setShape (pMesh m, pShape s, bool project=true);
int pumi_shape_getNumNode (pShape s, int topo);
bool pumi_shape_hasNode (pShape s, int topo);

void pumi_node_getCoord(pMeshEnt e, int i, double* xyz);
void pumi_node_setCoord(pMeshEnt e, int i, double* xyz);
void pumi_node_getCoordVector(pMeshEnt e, int i, Vector3& xyz);
void pumi_node_setCoordVector(pMeshEnt e, int i, Vector3 const& xyz);
/** double[3] cross product */
Vector3 pumi_vector3_cross(Vector3 const& a, Vector3 const& b);

pShape pumi_shape_getLagrange (int order);
pShape pumi_shape_getSerendipity ();
pShape pumi_shape_getConstant (int dimension);
pShape pumi_shape_getIP (int dimension, int order);
pShape pumi_shape_getVoronoi (int dimension, int order);
pShape pumi_shape_getIPFit(int dimension, int order);
pShape pumi_shape_getHierarchic (int order);

//************************************
//  Node numbering
//************************************
pGlobalNumbering pumi_numbering_createGlobal(pMesh m, const char* name, 
                 pShape shape=NULL, int num_component=1);
void pumi_numbering_deleteGlobal(pGlobalNumbering gn);
int pumi_mesh_getNumGlobalNumbering (pMesh m);
pGlobalNumbering pumi_mesh_getGlobalNumbering (pMesh m, int i);

pNumbering pumi_numbering_create (pMesh m, const char* name, pShape shape=NULL, int num_component=1);
pNumbering pumi_numbering_createLocalNode (pMesh m, const char* name, pShape shape=NULL);
pNumbering pumi_numbering_createOwned (pMesh m, const char* name, int dim);
pNumbering pumi_numbering_createOwnedNode (pMesh m, const char* name, pShape shape=NULL);
void pumi_numbering_delete(pNumbering n);
int pumi_numbering_getNumNode(pNumbering n);

//************************************
//  Field Management
//************************************

pField pumi_field_create(pMesh m, const char* name,
    int num_dof_per_node, int field_type=PUMI_PACKED, pShape shape = NULL);
int pumi_field_getSize(pField f);
int pumi_field_getType(pField f);
std::string pumi_field_getName(pField f);
pShape pumi_field_getShape (pField f);
void pumi_field_delete(pField f);
void pumi_field_synchronize(pField f);
void pumi_field_accumulate(pField f);
void pumi_field_freeze(pField f);
void pumi_field_unfreeze(pField f);
pField pumi_mesh_findField(pMesh m, const char* name);
int pumi_mesh_getNumField(pMesh m);
pField pumi_mesh_getField(pMesh m, int i);
void pumi_ment_getField (pMeshEnt e, pField f, int i, double* dof_data);
void pumi_ment_setField (pMeshEnt e, pField f, int i, double* dof_data);
// verify field
void pumi_field_verify(pMesh m, pField f=NULL);
void pumi_field_print(pField f);
#endif
