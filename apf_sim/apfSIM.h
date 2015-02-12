#ifndef APFSIM_H
#define APFSIM_H

#include <apfMesh.h>
#include <PartitionedMeshTypes.h>

namespace apf {
/** \brief Creates an apf::Mesh from a Simmetrix mesh.
  *
  * \details This object should be destroyed by apf::destroyMesh.
  */
Mesh* createMesh(pParMesh mesh);

/** \brief Casts a Simmetrix entity to an apf::MeshEntity.
  *
  * \details This does not create any objects, use freely.
  */
MeshEntity* castEntity(pEntity entity);

class TagSIM;

class MeshSIM : public Mesh
{
  public:
    MeshSIM(pParMesh m);
    virtual ~MeshSIM();
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
    int getType(MeshEntity* e);
    void getRemotes(MeshEntity* e, Copies& remotes);
    void getResidence(MeshEntity* e, Parts& residence);
    MeshTag* createDoubleTag(const char* name, int size);
    MeshTag* createIntTag(const char* name, int size);
    MeshTag* createLongTag(const char* name, int size);
    MeshTag* findTag(const char* name);
    void destroyTag(MeshTag* tag);
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

}//namespace apf

#endif
