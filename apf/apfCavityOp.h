/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFCAVITYOP_H
#define APFCAVITYOP_H

/** \file apfCavityOp.h
  \brief Parallel mesh cavity operator system */

#include "apfMesh.h"
#include <vector>
#include <cstring>

namespace apf {

/** \page cavity Cavity Operator
   A cavity operation is a generalization for
   anything operating on a mesh entity
   that needs access to a surrounding cavity
   comprised of adjacent elements and their
   closures, either to read or write attached
   data or even to modify the mesh.

   A user-provided apf::CavityOp should define two
   methods: setEntity() and apply().
   setEntity() should proceed as follows:
     1) check if the operator should be applied.
        if not, return SKIP.
     2) try to build the cavity. before getting
        upward adjacencies of entities, call
        CavityOp::requestLocality on them. If it returns
        false, return REQUEST.
     3) return OK

   Entity x needs to be local before upward
   adjacent entities of x are added to the cavity.
   If setEntity returns OK, APF may subsequently
   call apply() which should do all the actual
   writing or mesh modification.
   setEntity should not write or modify the cavity.
   if setEntity returns REQUEST, it may
   be called several times.
   Any call to requestLocality which returns
   false will request a migration to make all the
   given entities local.
   It is for correctness to request all the
   surface entities at some stage of cavity expansion
   before returning REQUEST.

   To use a CavityOp, currently the function
   applyToDimension is provided which will
   apply the operator where needed to all mesh
   entities of a given dimension in the entire
   distributed mesh.

   To have an efficient CavityOp, setEntity() should
   store the cavity as a local variable for apply() to use.

   mesh modifying operators should call preDeletion(e) before
   actually deleting an entity to prevent a crash due to
   iterator invalidation.
*/

/** \brief user-defined mesh cavity operator */
class CavityOp
{
  public:
    /** \brief constructor
      \param canModify true iff the operator can create or
                       destroy mesh entities */
    CavityOp(Mesh* m, bool canModify = false);
    /** \brief outcome of a setEntity call */
    enum Outcome {
      /** \brief skip the given entity */
      SKIP,
      /** \brief apply to the given entity */
      OK,
      /** \brief request more elements around the entity */
      REQUEST };
    /** \brief evaluate whether what to do with this entity */
    virtual Outcome setEntity(MeshEntity* e) = 0;
    /** \brief apply the operator on the (now local) cavity */
    virtual void apply() = 0;
    /** \brief parallel collective operation over entities of one dimension */
    void applyToDimension(int d);
    /** \brief within setEntity, require that entities be made local */
    bool requestLocality(MeshEntity** entities, int count);
    /** \brief call before deleting a mesh entity during the operation */
    void preDeletion(MeshEntity* e);
    /** \brief mesh pointer for convenience */
    Mesh* mesh;
  private:
    typedef std::vector<MeshEntity*> Requests;
    Requests requests;
    bool isRequesting;
    struct PullRequest { MeshEntity* e; int to; };
    bool sendPullRequests(std::vector<PullRequest>& received);
    bool tryToPull();
    void applyLocallyWithModification(int d);
    void applyLocallyWithoutModification(int d);
    bool canModify;
    bool movedByDeletion;
    MeshIterator* iterator;
};

} //namespace apf

#endif
