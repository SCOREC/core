/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFMESH_H
#define APFMESH_H

#include <vector>
#include <map>
#include <set>
#include "apfVector.h"
#include "apfDynamicArray.h"

namespace apf {

class FieldShape;
class Field;
template <class T>
class NumberingOf;
typedef NumberingOf<int> Numbering;

class MeshEntity;
class MeshIterator;
class MeshTag;
class ModelEntity;

typedef std::map<int,MeshEntity*> Copies;
typedef std::set<int> Parts;
typedef DynamicArray<MeshEntity*> Adjacent;
/* a static array type for higher-performance downward adjacency queries */
typedef MeshEntity* Downward[12];

class Migration;

struct Up
{
  int n;
  MeshEntity* e[256];
};

struct Match
{
  int peer;
  MeshEntity* entity;
};

typedef DynamicArray<Match> Matches;

class Mesh
{
  public:
    void init(FieldShape* s);
    virtual ~Mesh();
    /* returns the element dimension of this mesh */
    virtual int getDimension() = 0;
    /* returns the number of entities in this dimension */
    virtual std::size_t count(int dimension) = 0;
    /* begins iteration over elements of one dimension */
    virtual MeshIterator* begin(int dimension) = 0;
    /* iterate over mesh entities
       0 is returned at the end of the iteration */
    virtual MeshEntity* iterate(MeshIterator* it) = 0;
    /* destroy the iterator. and end() call should match
       every begin() call to prevent memory leaks */
    virtual void end(MeshIterator* it) = 0;
    /* Returns true if the entity is shared in parallel */
    virtual bool isShared(MeshEntity* e) = 0;
    /* Returns true if the entity is shared in parallel
       and this is the dominant copy, or the entity is not shared. */
    virtual bool isOwned(MeshEntity* e) = 0;
    /* Returns the owning part number of this entity */
    virtual int getOwner(MeshEntity* e) = 0;
    enum { VERTEX, EDGE, TRIANGLE, QUAD, TET, HEX, PRISM, PYRAMID, TYPES };
    static int const adjacentCount[TYPES][4];
    static int const typeDimension[TYPES];
    /* Returns the set of entities of one dimension adjacent to a
       given entity. Unlike some databases, this includes the
       entity itself if the same dimension is given. */
    virtual void getAdjacent(MeshEntity* e, int dimension, Adjacent& adjacent) = 0;
    virtual int getDownward(MeshEntity* e, int dimension, MeshEntity** adjacent) = 0;
    virtual int countUpward(MeshEntity* e) = 0;
    virtual MeshEntity* getUpward(MeshEntity* e, int i) = 0;
    virtual void getUp(MeshEntity* e, Up& up) = 0;
    virtual bool hasUp(MeshEntity* e) = 0;
    /* Returns the coordinates of a mesh coordinate field node
       on a mesh entity. Most databases support at most
       one node on a vertex and one node on an edge for
       a 2nd-order Serendipity coordinate field. */
    void getPoint(MeshEntity* e, int node, Vector3& point);
    virtual void getPoint_(MeshEntity* e, int node, Vector3& point) = 0;
    virtual void getParam(MeshEntity* e, Vector3& p) = 0;
    /* Returns the type of a mesh entity using the Mesh::Type enumeration */
    virtual int getType(MeshEntity* e) = 0;
    virtual void getRemotes(MeshEntity* e, Copies& remotes) = 0;
    virtual void getResidence(MeshEntity* e, Parts& residence) = 0;
    /* Creates a double array tag over the mesh given a name and size */
    virtual MeshTag* createDoubleTag(const char* name, int size) = 0;
    /* Creates an int array tag over the mesh given a name and size */
    virtual MeshTag* createIntTag(const char* name, int size) = 0;
    /* Creates a long array tag over the mesh given a name and size */
    virtual MeshTag* createLongTag(const char* name, int size) = 0;
    /* Finds a tag by name, returns 0 if it doesn't exist */
    virtual MeshTag* findTag(const char* name) = 0;
    /* Removes a mesh tag. This does not detach data from entities. */
    virtual void destroyTag(MeshTag* tag) = 0;
    virtual void getTags(DynamicArray<MeshTag*>& tags) = 0;
    /* get/set double array tag data */
    virtual void getDoubleTag(MeshEntity* e, MeshTag* tag, double* data) = 0;
    virtual void setDoubleTag(MeshEntity* e, MeshTag* tag, double const* data) = 0;
    /* get/set int array tag data */
    virtual void getIntTag(MeshEntity* e, MeshTag* tag, int* data) = 0;
    virtual void setIntTag(MeshEntity* e, MeshTag* tag, int const* data) = 0;
    /* get/set long array tag data */
    virtual void getLongTag(MeshEntity* e, MeshTag* tag, long* data) = 0;
    virtual void setLongTag(MeshEntity* e, MeshTag* tag, long const* data) = 0;
    /* Detach tag data from an entity. must be attached already */
    virtual void removeTag(MeshEntity* e, MeshTag* tag) = 0;
    /* Returns true if there is data for this tag attached */
    virtual bool hasTag(MeshEntity* e, MeshTag* tag) = 0;
    enum { DOUBLE, INT, LONG };
    virtual int getTagType(MeshTag* t) = 0;
    virtual int getTagSize(MeshTag* t) = 0;
    virtual const char* getTagName(MeshTag* t) = 0;
    /* return the model entity dimension */
    virtual int getModelType(ModelEntity* e) = 0;
    /* get the dimension-unique model entity identifier */
    virtual int getModelTag(ModelEntity* e) = 0;
    /* get the model entity by dimension and identifier */
    virtual ModelEntity* findModelEntity(int type, int tag) = 0;
    /* get geometric classification */
    virtual ModelEntity* toModel(MeshEntity* e) = 0;
    virtual void getModelFaceNormal(ModelEntity* face, Vector3 const& p, Vector3& n) = 0;
    virtual void getModelEdgeTangent(ModelEntity* edge, double p, Vector3& n) = 0;
    /* evaluate parametric coordinate (p) as a spatial point (x)
       Return true iff the evaluation worked (i.e. the modeler supports it) */
    virtual bool snapToModel(ModelEntity* m, Vector3 const& p, Vector3& x) = 0;
    /* static lookup table for the dimension of each type in
       Mesh::Type */
    static int getEntityDimension(int type);
    /* get the shape of the mesh's coordinate field */
    FieldShape* getShape() const;
    Field* getCoordinateField() {return coordinateField;}
    void changeCoordinateField(Field* f);
    /* Migrate elements.
       This is called with a mapping from elements
       to process IDs, which will be deleted during migration */
    virtual void migrate(Migration* plan) = 0;
    virtual int getId() = 0;
    /* write the underlying mesh into a set of files */
    virtual void writeNative(const char* fileName) = 0;
    /* actually destroy the underlying mesh data structure */
    virtual void destroyNative() = 0;
    /* run a set of consistency checks on the underlying data structure */
    virtual void verify() = 0;
    virtual bool hasMatching() = 0;
    virtual void getMatches(MeshEntity* e, Matches& m) = 0;
    void addField(Field* f);
    void removeField(Field* f);
    Field* findField(const char* name);
    int countFields();
    Field* getField(int i);
    void addNumbering(Numbering* f);
    void removeNumbering(Numbering* f);
    Numbering* findNumbering(const char* name);
    int countNumberings();
    Numbering* getNumbering(int i);
    bool hasFrozenFields;
  protected:
    Field* coordinateField;
    std::vector<Field*> fields;
    std::vector<Numbering*> numberings;
};

extern int const tri_edge_verts[3][2];
extern int const quad_edge_verts[4][2];
extern int const tet_edge_verts[6][2];
extern int const tet_tri_verts[4][3];
extern int const prism_tri_verts[2][3];
extern int const prism_quad_verts[3][4];
extern int const pyramid_tri_verts[4][3];

void unite(Parts& into, Parts const& from);

void getFacePeers(Mesh* m, Parts& peers);
/* given a mesh iterator from Mesh::begin, will use it to iterate
   only over entities of that dimension which have a copy on
   the given part. Call Mesh::end after as usual. */
MeshEntity* iterateBoundary(Mesh* m, MeshIterator* it, int part);

class Migration
{
  public:
    Migration(Mesh* m);
    ~Migration();
    int count();
    MeshEntity* get(int i);
    bool has(MeshEntity* e);
    void send(MeshEntity* e, int to);
    int sending(MeshEntity* e);
    Mesh* getMesh() {return mesh;}
  private:
    Mesh* mesh;
    MeshTag* tag;
    std::vector<MeshEntity*> elements;
};

void removeTagFromDimension(Mesh* m, MeshTag* tag, int d);

int findIn(MeshEntity** a, int n, MeshEntity* e);

/* finds the entity with the given downward adjacencies,
   regardless of ordering */
MeshEntity* findUpward(Mesh* m, int type, MeshEntity** down);

void findTriDown(
    Mesh* m,
    MeshEntity** verts,
    MeshEntity** down);
/* finds an element from a set of vertices */
MeshEntity* findElement(
    Mesh* m,
    int type,
    MeshEntity** verts);

MeshEntity* getEdgeVertOppositeVert(Mesh* m, MeshEntity* edge, MeshEntity* v);

int countEntitiesOfType(Mesh* m, int type);

class ElementVertOp
{
  public:
    virtual MeshEntity* apply(int type, MeshEntity** down) = 0;
    MeshEntity* run(int type, MeshEntity** verts);
    void runDown(int type, MeshEntity** verts, MeshEntity** down);
};

void changeMeshShape(Mesh* m, FieldShape* newShape, bool project = true);

void unfreezeFields(Mesh* m);

std::pair<int,MeshEntity*> getOtherCopy(Mesh* m, MeshEntity* s);

/* equivalent to m->isShared(e) if there is no matching.
   Otherwise, returns whether the entity has remote copies or matched copies */
bool hasCopies(Mesh* m, MeshEntity* e);
/* equivalent to m->isOwned(e) if there is no matching.
   Otherwise, implements an ownership scheme over matched entities */
bool isOriginal(Mesh* m, MeshEntity* e);

int getDimension(Mesh* m, MeshEntity* e);

void verify(Mesh* m);

bool isSimplex(int type);

} //namespace apf

#endif
