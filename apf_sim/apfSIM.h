#ifndef APFSIM_H
#define APFSIM_H

#include <MeshSim.h>

#include <apfMesh2.h>
#include <PartitionedMeshTypes.h>

void P_setPos(pPoint,double,double,double);
class Field;

namespace apf {

/**
 * \brief Creates an apf::Mesh from a Simmetrix mesh.
 *
 * \details This object should be destroyed by apf::destroyMesh.
 */
Mesh2* createMesh(pParMesh mesh);

/**
  * \brief Casts a Simmetrix entity to an apf::MeshEntity.
  *
  * \details This does not create any objects, use freely.
  */
MeshEntity* castEntity(pEntity entity);

class TagSIM;

class MeshSIM : public Mesh2
{
  public:
    MeshSIM(pParMesh m);
    virtual ~MeshSIM();
    // Mesh interface =======================
    int getDimension();
    std::size_t count(int dimension);
    MeshIterator* begin(int dimension);
    MeshEntity* iterate(MeshIterator* it);
    void end(MeshIterator* it);
    bool isShared(MeshEntity* e);
    bool isOwned(MeshEntity* e);
    int getOwner(MeshEntity* e);
    void getAdjacent(MeshEntity* e, int dimension, DynamicArray<MeshEntity*>& adjacent);
    int getDownward(MeshEntity* e, int dimension, MeshEntity** adjacent);
    int countUpward(MeshEntity* e);
    MeshEntity* getUpward(MeshEntity* e, int i);
    virtual void getUp(MeshEntity* e, Up& up);
    virtual bool hasUp(MeshEntity* e);
    void getPoint_(MeshEntity* e, int node, Vector3& point);
    void getParam(MeshEntity* e, Vector3& point);
    Type getType(MeshEntity* e);
    void getRemotes(MeshEntity* e, Copies& remotes);
    void getResidence(MeshEntity* e, Parts& residence);
    MeshTag* createDoubleTag(const char* name, int size);
    MeshTag* createIntTag(const char* name, int size);
    MeshTag* createLongTag(const char* name, int size);
    MeshTag* findTag(const char* name);
    void destroyTag(MeshTag* tag);
    void renameTag(MeshTag*, const char*);
    unsigned getTagChecksum(MeshTag*,int);
    void getTags(DynamicArray<MeshTag*>& ts);
    void getDoubleTag(MeshEntity* e, MeshTag* tag, double* data);
    void setDoubleTag(MeshEntity* e, MeshTag* tag, double const* data);
    void getIntTag(MeshEntity* e, MeshTag* tag, int* data);
    void setIntTag(MeshEntity* e, MeshTag* tag, int const* data);
    void getLongTag(MeshEntity* e, MeshTag* tag, long* data);
    void setLongTag(MeshEntity* e, MeshTag* tag, long const* data);
    void removeTag(MeshEntity* e, MeshTag* tag);
    bool hasTag(MeshEntity* e, MeshTag* tag);
    int getTagType(MeshTag* t);
    int getTagSize(MeshTag* t);
    const char* getTagName(MeshTag* t);
    ModelEntity* toModel(MeshEntity* e);
    gmi_model* getModel();
    void migrate(Migration* plan);
    int getId();
    void writeNative(const char* fileName);
    void destroyNative();
    void verify();
    bool hasMatching() {return hasMatches;}
    void getMatches(MeshEntity* e, Matches& m);
    bool hasMatches;
    bool isGhosted(MeshEntity*) { return false; }
    bool isGhost(MeshEntity*) { return false; }
    int getGhosts(MeshEntity*, Copies&) { return 0; }
    // Mesh2 interface ==============================
    void setRemotes(MeshEntity*, Copies&) {}
    void addRemote(MeshEntity*, int, MeshEntity*) {}
    void setResidence(MeshEntity*, Parts&) {}
    void setParam(MeshEntity*, Vector3 const &) {}
    void increment(MeshIterator*) {}
    bool isDone(MeshIterator*) { return true; }
    MeshEntity * deref(MeshIterator*) { return NULL; }
    void setPoint_(MeshEntity* me, int node, Vector3 const & p);
    MeshEntity * createVert_(ModelEntity*) { return NULL; }
    MeshEntity * createEntity_(int, ModelEntity*, MeshEntity**) { return NULL; }
    void destroy_(MeshEntity* ) {}
    void setModelEntity(MeshEntity*, ModelEntity*) {}
    void addMatch(MeshEntity*, int, MeshEntity* ) {}
    void clearMatches(MeshEntity*) {}
    void clear_() {}
    void acceptChanges() {}
    void addGhost(MeshEntity*, int, MeshEntity*) {}
    void deleteGhost(MeshEntity*) {}
    pParMesh getMesh() { return mesh; }
  protected:
    pParMesh mesh;
    pMesh part;
    int iterDim;
    VIter viter;
    EIter eiter;
    FIter fiter;
    RIter riter;
    int d;
    std::vector<TagSIM*> tags;
    gmi_model* model;
};

apf::Field* createSIMGeneralField(Mesh* m, const char* name, int valueType, int components, FieldShape* shape);
::Field* getSIMField(apf::Field* f);
apf::Field* wrapSIMField(Mesh* m, ::Field* fd);
apf::Field* createSIMField(Mesh* m, const char* name, int valueType, FieldShape* shape);
apf::Field* createSIMLagrangeField(Mesh* m, const char* name, int valueType, int order);
apf::Field* createSIMFieldOn(Mesh* m, const char* name, int valueType);
apf::Field* createSIMPackedField(Mesh* m, const char* name, int components, apf::FieldShape* shape = 0);

template <typename T>
static void pListToDynamicArray(pPList list, DynamicArray<T>& array)
{
  int n = PList_size(list);
  array.setSize(n);
  for (int i=0; i < n; ++i)
    array[i] = reinterpret_cast<T>(PList_item(list,i));
  PList_delete(list);
}

template<typename T>
static void DynamicArrayTopList(const DynamicArray<T> & array, pPList & list)
{
  int n = array.getSize();
  PList_clear(list);
  for(int ii = 0; ii < n; ii++)
    list = PList_append(list,reinterpret_cast<void*>(array(ii)));
}



}//namespace apf

#endif
