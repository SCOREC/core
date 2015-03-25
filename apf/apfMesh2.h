/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_MESH2_H
#define APF_MESH2_H

/** \file apfMesh2.h
    \brief The APF Mesh modification interface */

#include "apfMesh.h"
#include <set>

namespace apf {

/** \brief Extended mesh interface for modification
  \details this interface, which is a superset of apf::Mesh,
  includes methods for mesh modification.
  Not all mesh databases support this, which is why the
  two classes are separated.

  In this case, mesh modification means any entity creation,
  deletion, changing inter-part boundary links (including
  remote copies, matching, etc), and changing node coordinates. */

class Mesh2 : public Mesh
{
  public:
/** \brief Set the remote copies of an entity
 \details this does not affect the residence or ownership,
 so users are advised to use apf::stitchMesh after calling
 this function. */
    virtual void setRemotes(MeshEntity* e, Copies& remotes) = 0;
/** \brief Add just one remote copy to an entity */
    virtual void addRemote(MeshEntity* e, int p, MeshEntity* r) = 0;
/** \brief Set the resident part set of an entity
  \details this is also known as partition model classification */
    virtual void setResidence(MeshEntity* e, Parts& residence) = 0;
/** \brief Set the geometric parametric coordinates for a vertex */
    virtual void setParam(MeshEntity* e, Vector3 const& p) = 0;
/** \brief Just increment an iterator */
    virtual void increment(MeshIterator* it) = 0;
/** \brief Return true iff the iterator points past the end */
    virtual bool isDone(MeshIterator* it) = 0;
/** \brief Just dereference an iterator without incrementing it
  \details this is needed by apf::CavityOp, and Simmetrix meshes
  don't support the separation of dereference and increment */
    virtual MeshEntity* deref(MeshIterator* it) = 0;
/** \brief Set the spacial coordinates of a mesh node */
    void setPoint(MeshEntity* e, int node, Vector3 const& p);
/** \brief Underlying implementation of apf::Mesh2::setPoint */
    virtual void setPoint_(MeshEntity* e, int node, Vector3 const& p) = 0;
/** \brief Create a fully-specified vertex
  \details this function sets geometric classification and coordinates
  all at once. see apf::Mesh2::createVert for a more minimal interface */
    MeshEntity* createVertex(ModelEntity* c, Vector3 const& point,
        Vector3 const& param);
/** \brief require that no fields are stored in arrays */
    void requireUnfrozen()
    {
      if (hasFrozenFields)
        unfreezeFields(this);
    }
/** \brief Underlying implementation of apf::Mesh2::createVert */
    virtual MeshEntity* createVert_(ModelEntity* c) = 0;
/** \brief Just create a vertex
  \param c geometric classification, which is very immutable in APF */
    MeshEntity* createVert(ModelEntity* c)
    {
      requireUnfrozen();
      return createVert_(c);
    }
/** \brief Underlying implementation of apf::Mesh2::createEntity */
    virtual MeshEntity* createEntity_(int type, ModelEntity* c,
        MeshEntity** down) = 0;
/** \brief Create a non-vertex mesh entity
  \details to create entities from more than one level down,
  including intermediate entities, see apf::buildElement
  \param type select from apf::Mesh::Type
  \param c geometric classification, which is very immutable in APF
  \param down array of one-level downward adjacent entities */
    MeshEntity* createEntity(int type, ModelEntity* c, MeshEntity** down)
    {
      requireUnfrozen();
      return createEntity_(type,c,down);
    }
/** \brief Underlying implementation of apf::Mesh2::destroy */
    virtual void destroy_(MeshEntity* e) = 0;
/** \brief Destroy a mesh entity
  \details this does not destroy any other entities, including
  downward adjacencies */
    void destroy(MeshEntity* e)
    {
      requireUnfrozen();
      destroy_(e);
    }
/** \brief Add a matched copy to an entity */
    virtual void addMatch(MeshEntity* e, int peer, MeshEntity* match) = 0;
/** \brief Remove all matched copies of an entity */
    virtual void clearMatches(MeshEntity* e) = 0;
/** \brief Implementation-defined synchronization after modification
  \details users are encouraged to call this function after finishing
  mesh modifications so that all structures are properly updated before
  using the mesh any further. */
    virtual void acceptChanges() = 0;
};

/** \brief APF's migration function, works on apf::Mesh2
 \details if your database implements apf::Mesh2
 (and residence is separate from remote copies)
 then you may use this to implement most of apf::Mesh::migrate.
 Users of APF are encouraged to call apf::Mesh::migrate instead
 of calling this directly */
void migrate(Mesh2* m, Migration* plan);

void migrateSilent(Mesh2* m, Migration* plan);

/** \brief set the maximum elements that apf::migrate moves at once
  \details apf::migrate implements gradual limited migration
  in an effort to help applications keep memory use to a minimum.
  This function globally sets the limit on migration, which
  causes any migration requests greater than the limit to
  be performed as several consecutive migrations. */
void setMigrationLimit(size_t maxElements);

class Field;

/** \brief add a field (times a factor) to the mesh coordinates
  \details this is useful in mechanical deformation problems
  to transform the mesh from reference space to deformed space.
  Setting the factor to -1 will undo the deformation */
void displaceMesh(Mesh2* m, Field* d, double factor=1.0);

/** \brief User-defined entity creation callback */
class BuildCallback
{
  public:
    /** \brief will be called after an entity is created */
    virtual void call(MeshEntity* e) = 0;
};

/** \brief like apf::Mesh2::createEntity,
  but returns already existing entities */
MeshEntity* makeOrFind(
    Mesh2* m,
    ModelEntity* c,
    int type,
    MeshEntity** down,
    BuildCallback* cb = 0);

/** \brief build an entity from its vertices
  \details any missing intermediate entities will also be built,
  and all will be classified to (c).
  If a non-zero BuildCallback pointer is given, it will be called
  for each entity created, including intermediate ones. */
MeshEntity* buildElement(
    Mesh2* m,
    ModelEntity* c,
    int type,
    MeshEntity** verts,
    BuildCallback* cb = 0);

/** \brief build a one-element mesh
  \details this is mostly useful for debugging
  \todo this doesn't get used much, maybe remove it */
MeshEntity* buildOneElement(
    Mesh2* m,
    ModelEntity* c,
    int type,
    Vector3 const* points);

/** \brief Set entity residence based on remote copies
  \details this function acts on all entities of one dimension */
void initResidence(Mesh2* m, int dim);

/** \brief infer all remote copies from those of vertices
  \details
  given that the remote copies of the vertices are set up correctly, this
  function will synchronize the remote copies and resident part sets for all
  other entities correctly. */
void stitchMesh(Mesh2* m);

void packDataClone(Mesh2* m, int to);
void unpackDataClone(Mesh2* m);

}//namespace apf

#endif
