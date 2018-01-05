/******************************************************************************

  Copyright 2013 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

*******************************************************************************/
#ifndef MA_FACE_SPLIT_H
#define MA_FACE_SPLIT_H

#include "maAdapt.h"
#include <pcu_util.h>

namespace ma {

/** \brief Single face split operator
 *  \detail Restricted to tetrahedral mesh only.
 */
class FaceSplit
{
  public:
    FaceSplit(Adapt* a);

    /** \brief Set the target face
     *  \param e target MeshEntity (face)
     *  \return whether or not the face was added
     */
    bool setFace(Entity* e);
    void makeNewElements();
    void cancel();

    /** \brief Perform solution transfer */
    void transfer();
    void destroyOldElements();

    /** Returns the vertex created in this split */
    Entity* getSplitVert();
    EntityArray& getTets() {return toSplit[3];}
    Adapt* getAdapt() {return adapter;}
  private:
    Adapt* adapter;
    bool shouldCollect[4];
    EntityArray toSplit[4];
    apf::DynamicArray<EntityArray> newEntities[4];
};

Entity* makeSplitVertOnFace(Adapt* a, Entity* face);

/** Split a triangular face with a vertex in the center */
Entity* splitTri0(Adapt* a, Entity* parent);

}
#endif
