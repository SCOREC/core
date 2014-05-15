/*
 * Copyright (C) 2011-2013 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFCOORDDDATA_H
#define APFCOORDDDATA_H

#include "apfFieldData.h"

namespace apf {

/* This is a special field storage implementation for
   the mesh's coordinate field.
   It pulls the data from apf::Mesh::getPoint_,
   which will take it from wherever the mesh database keeps
   its vertex and mid-edge node coordinates.

   So far it only supports linear and quadratic shapes since
   very high order shapes should really be handled on the APF
   side instead of copy-pasting APF functionality into a mesh
   library.
*/

class CoordData : public FieldDataOf<double>
{
  public:
    virtual void init(FieldBase* f);
    virtual bool hasEntity(MeshEntity* e);
    virtual void removeEntity(MeshEntity* e);
    virtual void get(MeshEntity* e, double* data);
    virtual void set(MeshEntity* e, double const* data);
    virtual bool isFrozen() { return false; }
  private:
    Mesh* mesh;
    FieldShape* shape;
};

}

#endif

