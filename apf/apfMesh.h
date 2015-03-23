/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_MESH_H
#define APF_MESH_H

/** \file apfMesh.h
    \brief The APF Mesh interface*/

#include <vector>
#include <map>
#include <set>
#include "apfVector.h"
#include "apfDynamicArray.h"

struct gmi_model;

namespace apf {

class FieldShape;
class Field;
template <class T>
class NumberingOf;
/** \brief Numbering is meant to be a 32-bit local numbering */
typedef NumberingOf<int> Numbering;
typedef NumberingOf<long> GlobalNumbering;

class MeshEntity;
class MeshIterator;
class MeshTag;

class ModelEntity;

/** \brief Remote copy container.
  \details the key is the part id, the value
  is the on-part pointer to the remote copy */
typedef std::map<int,MeshEntity*> Copies;
/** \brief Set of unique part ids */
typedef std::set<int> Parts;
/** \brief Set of adjacent mesh entities
  \details see also apf::Downward and apf::Up */
typedef DynamicArray<MeshEntity*> Adjacent;
/** \brief a static array type downward adjacency queries.
    \details using statically sized arrays saves time by
             avoiding dynamic allocation, and downward
             adjacencies have a guaranteed bound.
 */
typedef MeshEntity* Downward[12];

class Migration;

/** \brief statically sized container for upward adjacency queries.
    \details see apf::Downward for static size rationale.
    Although our algorithmic complexity proofs rely on upward
    adjacencies being bound by a constant, this constant has yet
    to be pinpointed. Some (bad) Simmetrix meshes have around 300
    edges per vertex, so we are now at 400.
  */
struct Up
{
  /** \brief actual number of adjacent entities */
  int n;
  /** \brief array containing pointers to adjacent entities */
  MeshEntity* e[400];
};

/** \brief a reference to an object representing the same entity
  \details This may represent a copy along a partition boundary
  or a matched copy for periodic meshes.
  */
struct Copy
{
  /** \brief required */
  Copy() {}
  /** \brief build from contents */
  Copy(int p, MeshEntity* e):peer(p),entity(e) {}
  /** \brief resident part of the copy object */
  int peer;
  /** \brief on-part pointer to the copy object */
  MeshEntity* entity;
};
/** \brief matches are just a special case of copies */
typedef Copy Match;

/** \brief a set of copies, possibly multiple copies per part */
typedef DynamicArray<Copy> CopyArray;
/** \brief a set of matched copies */
typedef CopyArray Matches;

/** \brief Interface to a mesh part
  \details This base class is the interface for almost all mesh
  operations in APF. Code that interacts with a mesh should do
  so through an apf::Mesh interface object.
  Mesh databases should derive an interface object and implement
  all pure virtual functions to be usable from APF.
 */
class Mesh
{
  public:
    /** \brief initialize the base class structures.
      \param s the field distribution of the coordinate field,
               apf::getLagrange(1) is a good default
      */
    void init(FieldShape* s);
    /** \brief destroy the base class structures.
        \details this does not destroy the underlying data
                 structure, use apf::Mesh::destroyNative for that.
      */
    virtual ~Mesh();
    /** \brief returns the element dimension of this mesh */
    virtual int getDimension() = 0;
    /** \brief returns the number of entities in this dimension */
    virtual std::size_t count(int dimension) = 0;
    /** \brief begins iteration over elements of one dimension */
    virtual MeshIterator* begin(int dimension) = 0;
    /** \brief iterate over mesh entities
        \details 0 is returned at the end of the iteration */
    virtual MeshEntity* iterate(MeshIterator* it) = 0;
    /** \brief destroy an iterator.
        \details an end() call should match every begin()
                 call to prevent memory leaks */
    virtual void end(MeshIterator* it) = 0;
    /** \brief Returns true if the entity is shared in parallel */
    virtual bool isShared(MeshEntity* e) = 0;
    /** \brief Returns true if the entity is shared in parallel
              and this is the dominant copy, or the entity is not shared. */
    virtual bool isOwned(MeshEntity* e) = 0;
    /** \brief Returns the owning part number of this entity */
    virtual int getOwner(MeshEntity* e) = 0;
    /** \brief Entity topological types */
    enum Type {
      /** \brief vertex */
      VERTEX,   //0
      /** \brief edge */
      EDGE,     //1
      /** \brief triangle */
      TRIANGLE, //2
      /** \brief quadrilateral (square) */
      QUAD,     //3
      /** \brief tetrahedron */
      TET,      //4
      /** \brief hexahedron (cube, brick) */
      HEX,      //5
      /** \brief triangular prism (wedge) */
      PRISM,    //6
      /** \brief quadrilateral pyramid */
      PYRAMID,  //7
      /** \brief placeholder to set array sizes */
      TYPES };  //8
    /** \brief for a given entity type,
               number of adjacent entities of a given dimension */
    static int const adjacentCount[TYPES][4];
    /** \brief for a given entity type, its dimension. */
    static int const typeDimension[TYPES];
    /** \brief name strings for apf::Mesh::Type */
    static char const* const typeName[TYPES];
    /** \brief Returns the set of entities of one dimension adjacent to a
       given entity.
       \details prefer to use apf::Mesh::getDownward and apf::Mesh::getUp
                when possible, this function is only superior for upward
                adjacencies of more than one level.
      */
    virtual void getAdjacent(MeshEntity* e, int dimension,
        Adjacent& adjacent) = 0;
    /** \brief Returns an ordered set of downward adjacent entities.
        \details Downward adjacent entities follow a strict ordering
        which was defined at entity creation time and should be consistent
        with the topological orderings shown below:
        \image html region_faces.jpg
        \image html region_edges.jpg
        \param adjacent the output array. can be a user-sized static
               array if the size is known, otherwise use apf::Downward */
    virtual int getDownward(MeshEntity* e, int dimension,
        MeshEntity** adjacent) = 0;
    /** \brief Return the number of one-level upward adjacent entities. */
    virtual int countUpward(MeshEntity* e) = 0;
    /** \brief Get the i'th one-level upward adjacent entity. */
    virtual MeshEntity* getUpward(MeshEntity* e, int i) = 0;
    /** \brief Get the unordered set of one-level upward entities. */
    virtual void getUp(MeshEntity* e, Up& up) = 0;
    /** \brief Return true iff the entity has upward adjacencies. */
    virtual bool hasUp(MeshEntity* e) = 0;
    /** \brief Returns the coordinates of a node on a mesh entity.
       \details Most databases support at most
       one node on a vertex and one node on an edge for
       a 2nd-order Serendipity coordinate field.
       \param e the entity associated with the node
       \param node in the case of curved meshes, which node on the entity
       \param point nodal coordinates
     */
    void getPoint(MeshEntity* e, int node, Vector3& point);
    /** \brief Implementation-defined code for apf::Mesh::getPoint */
    virtual void getPoint_(MeshEntity* e, int node, Vector3& point) = 0;
    /** \brief Get the geometric parametric coordinates of a vertex */
    virtual void getParam(MeshEntity* e, Vector3& p) = 0;
    /** \brief Get the topological type of a mesh entity.
      \returns a value from the apf::Mesh::Type enumeration */
    virtual int getType(MeshEntity* e) = 0;
    /** \brief Get the remote copies of an entity */
    virtual void getRemotes(MeshEntity* e, Copies& remotes) = 0;
    /** \brief Get the resident parts of an entity
      \details this includes parts with remote copies and the
               current part as well */
    virtual void getResidence(MeshEntity* e, Parts& residence) = 0;
    /** \brief Creates a double array tag over the mesh given a name and size */
    virtual MeshTag* createDoubleTag(const char* name, int size) = 0;
    /** \brief Creates an int array tag over the mesh given a name and size */
    virtual MeshTag* createIntTag(const char* name, int size) = 0;
    /** \brief Creates a long array tag over the mesh given a name and size */
    virtual MeshTag* createLongTag(const char* name, int size) = 0;
    /** \brief Finds a tag by name, returns 0 if it doesn't exist */
    virtual MeshTag* findTag(const char* name) = 0;
    /** \brief Removes a mesh tag. This does not detach data from entities. */
    virtual void destroyTag(MeshTag* tag) = 0;
    /** \brief Get all the tags on a mesh part. */
    virtual void getTags(DynamicArray<MeshTag*>& tags) = 0;
    /** \brief get double array tag data */
    virtual void getDoubleTag(MeshEntity* e, MeshTag* tag, double* data) = 0;
    /** \brief set double array tag data */
    virtual void setDoubleTag(MeshEntity* e, MeshTag* tag, double const* data) = 0;
    /** \brief get int array tag data */
    virtual void getIntTag(MeshEntity* e, MeshTag* tag, int* data) = 0;
    /** \brief set int array tag data */
    virtual void setIntTag(MeshEntity* e, MeshTag* tag, int const* data) = 0;
    /** \brief get long array tag data */
    virtual void getLongTag(MeshEntity* e, MeshTag* tag, long* data) = 0;
    /** \brief set long array tag data */
    virtual void setLongTag(MeshEntity* e, MeshTag* tag, long const* data) = 0;
    /** \brief detach tag data from an entity.
        \details data must be attached already. */
    virtual void removeTag(MeshEntity* e, MeshTag* tag) = 0;
    /** \brief Returns true if there is data for this tag attached */
    virtual bool hasTag(MeshEntity* e, MeshTag* tag) = 0;
    /** \brief Tag data type enumeration */
    enum TagType {
      /** \brief 64-bit IEE754 floating-point number */
      DOUBLE,
      /** \brief signed 32-bit integer */
      INT,
      /** \brief signed 64-bit integer */
      LONG };
    /** \brief get the data type of a tag
        \return a value in apf::Mesh::TagType */
    virtual int getTagType(MeshTag* t) = 0;
    /** \brief return the array size of a tag */
    virtual int getTagSize(MeshTag* t) = 0;
    /** \brief return the name of a tag
      \returns a pointer to an internal C string.
               do not free this pointer */
    virtual const char* getTagName(MeshTag* t) = 0;
    /** \brief get geometric classification */
    virtual ModelEntity* toModel(MeshEntity* e) = 0;
    /** \brief get a GMI interface to the geometric model */
    virtual gmi_model* getModel() = 0;
    /** \brief return the model entity dimension */
    int getModelType(ModelEntity* e);
    /** \brief get the dimension-unique model entity identifier */
    int getModelTag(ModelEntity* e);
    /** \brief get the model entity by dimension and identifier */
    ModelEntity* findModelEntity(int type, int tag);
    /** \brief return true if the geometric model supports snapping */
    bool canSnap();
    /** \brief evaluate parametric coordinate (p) as a spatial point (x) */
    void snapToModel(ModelEntity* m, Vector3 const& p, Vector3& x);
    /** \brief reparameterize mesh vertex (e) onto model entity (g) */
    void getParamOn(ModelEntity* g, MeshEntity* e, Vector3& p);
    /** \brief get the periodic properties of a model entity
      \param range if periodic, the parametric range
      \returns true if (g) is periodic along this axis */
    bool getPeriodicRange(ModelEntity* g, int axis,
        double range[2]);
    /** \brief get the distribution of the mesh's coordinate field */
    FieldShape* getShape() const;
    /** \brief get the mesh's coordinate field */
    Field* getCoordinateField() {return coordinateField;}
    /** \brief replace the mesh's coordinate field */
    void changeCoordinateField(Field* f);
    /** \brief Migrate elements.
       \param plan a mapping from local elements
                   to part IDs, which will be deleted during migration */
    virtual void migrate(Migration* plan) = 0;
    /** \brief Get the part ID */
    virtual int getId() = 0;
    /** \brief write the underlying mesh into a set of files */
    virtual void writeNative(const char* fileName) = 0;
    /** \brief actually destroy the underlying mesh data structure */
    virtual void destroyNative() = 0;
    /** \brief run a set of consistency checks on the underlying
               data structure */
    virtual void verify() = 0;
    /** \brief return true if the mesh has matched entities */
    virtual bool hasMatching() = 0;
    /** \brief get the matches of an entity */
    virtual void getMatches(MeshEntity* e, Matches& m) = 0;
    /** \brief estimate mesh entity memory usage.
      \details this is used by Parma_WeighByMemory
      \param type a value from apf::Mesh::Type
      \returns an estimate of how many bytes are needed
      to store an entity of (type) */
    virtual double getElementBytes(int) {return 1.0;}
    /** \brief associate a field with this mesh
      \details most users don't need this, functions in apf.h
               automatically call it */
    void addField(Field* f);
    /** \brief disassociate a field from this mesh
      \details most users don't need this, functions in apf.h
               automatically call it */
    void removeField(Field* f);
    /** \brief lookup a field by its unique name */
    Field* findField(const char* name);
    /** \brief get the number of associated fields */
    int countFields();
    /** \brief get the i'th associated field */
    Field* getField(int i);
    /** \brief associate a numbering with this mesh
      \details most users don't need this, functions in apfNumbering.h
               automatically call it */
    void addNumbering(Numbering* f);
    /** \brief disassociate a numbering from this mesh
      \details most users don't need this, functions in apfNumbering.h
               automatically call it */
    void removeNumbering(Numbering* f);
    /** \brief lookup a numbering by its unique name */
    Numbering* findNumbering(const char* name);
    /** \brief get the number of associated numberings */
    int countNumberings();
    /** \brief get the i'th associated numbering */
    Numbering* getNumbering(int i);
    void addGlobalNumbering(GlobalNumbering* f);
    void removeGlobalNumbering(GlobalNumbering* f);
    int countGlobalNumberings();
    GlobalNumbering* getGlobalNumbering(int i);
    /** \brief true if any associated fields use array storage */
    bool hasFrozenFields;
  protected:
    Field* coordinateField;
    std::vector<Field*> fields;
    std::vector<Numbering*> numberings;
    std::vector<GlobalNumbering*> globalNumberings;
};

/** \brief run consistency checks on an apf::Mesh structure
  \details this can be used to implement apf::Mesh::verify.
  Other implementations may define their own. */
void verify(Mesh* m);

long verifyVolumes(Mesh* m, bool printVolumes = true);

/** \brief get the dimension of a mesh entity */
int getDimension(Mesh* m, MeshEntity* e);

/** \brief unite two sets of unique part ids
  \param into becomes the union */
void unite(Parts& into, Parts const& from);

/** \brief removes a tag from all entities of dimension (d) */
void removeTagFromDimension(Mesh* m, MeshTag* tag, int d);

/** \brief find an entity from one-level downward adjacencies
   \details this function ignores the ordering of adjacent entities */
MeshEntity* findUpward(Mesh* m, int type, MeshEntity** down);

/** \brief finds an entity from a set of vertices */
MeshEntity* findElement(
    Mesh* m,
    int type,
    MeshEntity** verts);

/** \brief get the other vertex of an edge */
MeshEntity* getEdgeVertOppositeVert(Mesh* m, MeshEntity* edge, MeshEntity* v);

/** \brief get 2nd-order adjacent entities */
void getBridgeAdjacent(Mesh* m, MeshEntity* origin,
    int bridgeDimension, int targetDimension, Adjacent& result);

/** \brief count all on-part entities of one topological type */
int countEntitiesOfType(Mesh* m, int type);

/** \brief return true if the topological type is a simplex */
bool isSimplex(int type);

/** \brief get the average of the entity's vertex coordinates
  \details this also works if given just a vertex,
  so its a convenient way to get the centroid of any entity */
Vector3 getLinearCentroid(Mesh* m, MeshEntity* e);

/** \brief Migration plan object: local elements to destinations. */
class Migration
{
  public:
/** \brief must be constructed with a mesh
  \details use (new apf::Migration(mesh)) to make these objects */
    Migration(Mesh* m);
    Migration(Mesh* m, MeshTag* existingTag);
    ~Migration();
/** \brief return the number of elements with assigned destinations */
    int count();
/** \brief get the i'th element with an assigned destination */
    MeshEntity* get(int i);
/** \brief return true if the element has been assigned a destination */
    bool has(MeshEntity* e);
/** \brief assign a destination part id to an element */
    void send(MeshEntity* e, int to);
/** \brief return the destination part id of an element */
    int sending(MeshEntity* e);
    Mesh* getMesh() {return mesh;}
  private:
    Mesh* mesh;
    MeshTag* tag;
    std::vector<MeshEntity*> elements;
};

/** \brief abstract description of entity copy sharing
  \details this interface abstracts over remote copies,
  matching, and possible user-defined sharing models.
  For example, users can define a new Sharing object
  that uses a different ownership rule. */
struct Sharing
{
  virtual ~Sharing() {}
/** \brief return true if the entity is owned */
  virtual bool isOwned(MeshEntity* e) = 0;
/** \brief get the copies of the entity */
  virtual void getCopies(MeshEntity* e,
      CopyArray& copies) = 0;
};

/** \brief create a default sharing object for this mesh
  \details for normal meshes, the sharing object just
  describes remote copies. For matched meshes, the
  sharing object describes matches for matched entities
  and remote copies for other entities */
Sharing* getSharing(Mesh* m);

/** \brief map from triangle edge order to triangle vertex order */
extern int const tri_edge_verts[3][2];
/** \brief map from quad edge order to quad vertex order */
extern int const quad_edge_verts[4][2];
/** \brief map from tet edge order to tet vertex order */
extern int const tet_edge_verts[6][2];
/** \brief map from prism edge order to prism vertex order */
extern int const prism_edge_verts[9][2];
/** \brief map from pyramid edge order to pyramid vertex order */
extern int const pyramid_edge_verts[8][2];
/** \brief map from tet triangle order to tet vertex order */
extern int const tet_tri_verts[4][3];
/** \brief map from prism triangle order to prism vertex order */
extern int const prism_tri_verts[2][3];
/** \brief map from prism quad order to prism vertex order */
extern int const prism_quad_verts[3][4];
/** \brief map from pyramid triangle order to pyramid vertex order */
extern int const pyramid_tri_verts[4][3];

/** \brief scan the part for [vtx|edge|face]-adjacent part ids */
void getPeers(Mesh* m, int d, Parts& peers);

/** \brief find pointer (e) in array (a) of length (n)
  \returns -1 if not found, otherwise i such that a[i] = e */
int findIn(MeshEntity** a, int n, MeshEntity* e);

/** \brief given the vertices of a triangle, find its edges
  \param down the resulting array of edges */
void findTriDown(
    Mesh* m,
    MeshEntity** verts,
    MeshEntity** down);

/** \brief change the distribution of the mesh coordinate nodes
    \param project whether to project coordinate values from the old field */
void changeMeshShape(Mesh* m, FieldShape* newShape, bool project = true);

/** \brief unfreeze all associated fields
  \details see apf::unfreezeField */
void unfreezeFields(Mesh* m);

/** \brief count the number of mesh entities classified on a model entity */
int countEntitiesOn(Mesh* m, ModelEntity* me, int dim);

/** \brief count the number of owned entities of dimension (dim) */
int countOwned(Mesh* m, int dim);

/** \brief print global mesh entity counts per dimension */
void printStats(Mesh* m);

/** \brief print to stderr the number of empty parts, if any */
void warnAboutEmptyParts(Mesh* m);

/** \brief given a mesh face, return its remote copy */
Copy getOtherCopy(Mesh* m, MeshEntity* s);

class ElementVertOp
{
  public:
    virtual MeshEntity* apply(int type, MeshEntity** down) = 0;
    MeshEntity* run(int type, MeshEntity** verts);
    void runDown(int type, MeshEntity** verts, MeshEntity** down);
};

/** \brief get the type of the first entity in this dimension */
int getFirstType(Mesh* m, int dim);

/** \brief boundary entity alignment to an element
  \param m the mesh
  \param elem the element
  \param boundary an entity on the boundary of elem
  \param which index of (boundary) in getDownward((elem)...)
  \param flip true iff orientation of (boundary) is opposite canonical
  \param rotate position of canonical vertex 0 in boundary vertices,
                or in boundary vertices reversed if (flip)==true */
void getAlignment(Mesh* m, MeshEntity* elem, MeshEntity* boundary,
    int& which, bool& flip, int& rotate);

void packString(std::string s, int to);
std::string unpackString();
void packTagInfo(Mesh* m, MeshTag* t, int to);
void unpackTagInfo(std::string& name, int& type, int& size);

} //namespace apf

#endif
