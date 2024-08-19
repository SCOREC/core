#ifndef APFCAP
#define APFCAP

#include <apfMesh2.h>
#include <PCU.h>

#include "CapstoneModule.h"
#include "CreateMG_Framework_Geometry.h"
#include "CreateMG_Framework_Mesh.h"

using namespace CreateMG;
using namespace CreateMG::Attribution;
using namespace CreateMG::Mesh;
using namespace CreateMG::Geometry;



class capEntity;
class capMesh;

namespace apf {

/**
 * \brief Creates an apf::Mesh from a CapStone mesh.
 *
 * \details This object should be destroyed by apf::destroyMesh.
 */
Mesh2* createMesh(MeshDatabaseInterface* mdb, GeometryDatabaseInterface* gdb, pcu::PCU *PCUObj);

/**
  * \brief Casts a CapStone entity to an apf::MeshEntity.
  *
  * \details This does not create any objects, use freely.
  */
MeshEntity* castEntity(capEntity* entity);

class TagCAP;

class MeshCAP : public Mesh2
{
  public:
    MeshCAP(MeshDatabaseInterface* mdb, GeometryDatabaseInterface* gdb);
    virtual ~MeshCAP();
    /* --------------------------------------------------------------------- */
    /* Category 00: General Mesh APIs */
    /* --------------------------------------------------------------------- */
    // REQUIRED Member Functions //
    int getDimension();
    std::size_t count(int dimension);
    Type getType(MeshEntity* e);
    void verify();
    // OPTIONAL Member Functions //
    void writeNative(const char* fileName);
    void destroyNative();

    /* --------------------------------------------------------------------- */
    /* Category 01: General Getters and Setters for vertex coordinates */
    /* --------------------------------------------------------------------- */
    // REQUIRED Member Functions //
    void getPoint_(MeshEntity* e, int, Vector3& point);
    void setPoint_(MeshEntity* e, int, Vector3 const& p);
    void getParam(MeshEntity* e, Vector3& p);
    void setParam(MeshEntity* e, Vector3 const& point);

    /* --------------------------------------------------------------------- */
    /* Category 02: Iterators */
    /* --------------------------------------------------------------------- */
    // REQUIRED Member Functions //
    MeshIterator* begin(int dimension);
    MeshEntity* iterate(MeshIterator* it);
    void end(MeshIterator* it);
    void increment(MeshIterator* it);
    bool isDone(MeshIterator* it);
    MeshEntity* deref(MeshIterator* it);

    /* --------------------------------------------------------------------- */
    /* Category 03: Adjacencies */
    /* --------------------------------------------------------------------- */
    // REQUIRED Member Functions //
    void getAdjacent(MeshEntity* e, int dimension, Adjacent& adjacent);
    int getDownward(MeshEntity* e, int dimension, MeshEntity** adjacent);
    MeshEntity* getUpward(MeshEntity* e, int i);
    bool hasUp(MeshEntity* e);
    // OPTIONAL Member Functions //
    bool hasAdjacency(int from_dim, int to_dim);
    void createAdjacency(int from_dim, int to_dim);
    void deleteAdjacency(int from_dim, int to_dim);
    void getUp(MeshEntity* e, Up& up);
    int countUpward(MeshEntity* e);

    /* --------------------------------------------------------------------- */
    /* Category 04: CAD model inquires */
    /* --------------------------------------------------------------------- */
    // REQUIRED Member Functions //
    ModelEntity* toModel(MeshEntity* e);
    // OPTIONAL Member Functions //
    gmi_model* getModel();
    void setModel(gmi_model* newModel);
    void setModelEntity(MeshEntity* e, ModelEntity* me);

    /* --------------------------------------------------------------------- */
    /* Category 05: Entity Creation/Deletion */
    /* --------------------------------------------------------------------- */
    // REQUIRED Member Functions //
    MeshEntity* createVert_(ModelEntity* me);
    MeshEntity* createEntity_(int type, ModelEntity* me, MeshEntity** down);
    void destroy_(MeshEntity* e);

    /* --------------------------------------------------------------------- */
    /* Category 06: Attachable Data Functionality */
    /* --------------------------------------------------------------------- */
    // REQUIRED Member Functions //
    MeshTag* createDoubleTag(const char* name, int size);
    MeshTag* createIntTag(const char* name, int size);
    MeshTag* createLongTag(const char* name, int size);
    MeshTag* findTag(const char* name);
    void destroyTag(MeshTag* t);
    void renameTag(MeshTag* t, const char* name);
    void getTags(DynamicArray<MeshTag*>& tags);
    /* void getTag(MeshEntity* e, MeshTag* t, void* data); */
    /* void setTag(MeshEntity* e, MeshTag* t, void const* data); */
    void getDoubleTag(MeshEntity* e, MeshTag* tag, double* data);
    void setDoubleTag(MeshEntity* e, MeshTag* tag, double const* data);
    void getIntTag(MeshEntity* e, MeshTag* tag, int* data);
    void setIntTag(MeshEntity* e, MeshTag* tag, int const* data);
    void getLongTag(MeshEntity* e, MeshTag* tag, long* data);
    void setLongTag(MeshEntity* e, MeshTag* tag, long const* data);
    void removeTag(MeshEntity* e, MeshTag* t);
    bool hasTag(MeshEntity* e, MeshTag* t);
    int getTagType(MeshTag* t);
    int getTagSize(MeshTag* t);
    const char* getTagName(MeshTag* t);
    unsigned getTagChecksum(MeshTag*,int);


    /* --------------------------------------------------------------------- */
    /* Category 07: Distributed Meshes */
    /* --------------------------------------------------------------------- */
    // REQUIRED Member Functions //
    bool isShared(MeshEntity* e);
    bool isGhost(MeshEntity*) { return false; }
    bool isGhosted(MeshEntity*) { return false; }
    bool isOwned(MeshEntity* e);
    int getOwner(MeshEntity* e);
    void getRemotes(MeshEntity* e, Copies& remotes);
    void getResidence(MeshEntity* e, Parts& residence);
    int getId();
    void setResidence(MeshEntity*, Parts&) {}
    void acceptChanges() {}
    // OPTIONAL Member Functions //
    void deleteGhost(MeshEntity*) {}
    void addGhost(MeshEntity*, int, MeshEntity*) {}
    int getGhosts(MeshEntity*, Copies&) { return 0; }
    void migrate(Migration* plan);
    void setRemotes(MeshEntity*, Copies&) {}
    void addRemote(MeshEntity*, int, MeshEntity*) {}
    void clearRemotes(MeshEntity*) {}


    /* --------------------------------------------------------------------- */
    /* Category 08: Periodic Meshes */
    /* --------------------------------------------------------------------- */
    // REQUIRED Member Functions //
    bool hasMatching() { return false; }
    void getMatches(MeshEntity* e, Matches& m);
    // OPTIONAL Member Functions //
    void addMatch(MeshEntity*, int, MeshEntity*) {}
    void clearMatches(MeshEntity*) {}
    void clear_() {}
    void getDgCopies(MeshEntity* e, DgCopies& dgCopies, ModelEntity* me);



    MeshDatabaseInterface* getMesh() { return meshInterface; }
  protected:
    /* CapstoneModule capModule; */
    MeshDatabaseInterface* meshInterface;
    GeometryDatabaseInterface  *geomInterface;
    /* AppContext                 *c; */
    int iterDim;
    int d;
    gmi_model* model;
    std::vector<TagCAP*> tags;
};


}//namespace apf

#endif
