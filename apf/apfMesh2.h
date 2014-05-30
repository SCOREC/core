/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFMESH2_H
#define APFMESH2_H

#include "apfMesh.h"
#include <set>

namespace apf {

class Mesh2 : public Mesh
{
  public:
    virtual void setRemotes(MeshEntity* e, Copies& remotes) = 0;
    virtual void addRemote(MeshEntity* e, int p, MeshEntity* r) = 0;
    virtual void setResidence(MeshEntity* e, Parts& residence) = 0;
    virtual void increment(MeshIterator* it) = 0;
    virtual bool isDone(MeshIterator* it) = 0;
    virtual MeshEntity* deref(MeshIterator* it) = 0;
    virtual bool canSnap() = 0;
    virtual void setParam(MeshEntity* e, Vector3 const& p) = 0;
    virtual void getParamOn(ModelEntity* g, MeshEntity* e, Vector3& p) = 0;
    virtual bool getPeriodicRange(ModelEntity* g, int axis, double range[2]) = 0;
    void setPoint(MeshEntity* e, int node, Vector3 const& p);
    virtual void setPoint_(MeshEntity* e, int node, Vector3 const& p) = 0;
    MeshEntity* createVertex(ModelEntity* c, Vector3 const& point, Vector3 const& param);
    void requireUnfrozen()
    {
      if (hasFrozenFields)
        unfreezeFields(this);
    }
    virtual MeshEntity* createVert_(ModelEntity* c) = 0;
    MeshEntity* createVert(ModelEntity* c)
    {
      requireUnfrozen();
      return createVert_(c);
    }
    virtual MeshEntity* createEntity_(int type, ModelEntity* c, MeshEntity** down) = 0;
    MeshEntity* createEntity(int type, ModelEntity* c, MeshEntity** down)
    {
      requireUnfrozen();
      return createEntity_(type,c,down);
    }
    virtual void destroy_(MeshEntity* e) = 0;
    void destroy(MeshEntity* e)
    {
      requireUnfrozen();
      destroy_(e);
    }
    virtual void addMatch(MeshEntity* e, int peer, MeshEntity* match) = 0;
    virtual void clearMatches(MeshEntity* e) = 0;
    virtual void repartition(MeshTag* elementWeights, double maximumImbalance) = 0;
    virtual void acceptChanges() = 0;
};

/* apf's custom migration function, will work for anyone
   implementing Mesh2 */
void migrate(Mesh2* m, Migration* plan);
void setMigrationLimit(size_t maxElements);

class Field;

void displaceMesh(Mesh2* m, Field* d, double factor=1.0);

class BuildCallback
{
  public:
    virtual void call(MeshEntity* e) = 0;
};

MeshEntity* makeOrFind(
    Mesh2* m,
    ModelEntity* c,
    int type,
    MeshEntity** down,
    BuildCallback* cb = 0);

/* constructs a new element, any entities created
   will be classified to "c" */
MeshEntity* buildElement(
    Mesh2* m,
    ModelEntity* c,
    int type,
    MeshEntity** verts,
    BuildCallback* cb = 0);

MeshEntity* buildOneElement(
    Mesh2* m,
    ModelEntity* c,
    int type,
    Vector3* points);

void initResidence(Mesh2* m, int dim);
/* given that the remote copies of the vertices are set up correctly, this
   function will synchronize the remote copies and resident part sets for all
   other entities correctly */
void stitchMesh(Mesh2* m);

}//namespace apf

#endif
