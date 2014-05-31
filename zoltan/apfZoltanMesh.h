/*
 * Copyright (C) 2014 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_ZOLTAN_MESH_H
#define APF_ZOLTAN_MESH_H

#include <apfMesh.h>
#include <apfNumbering.h>

namespace apf {

class ZoltanMesh
{
  public:
    ZoltanMesh(Mesh* mesh_, bool local, int method_, int approach_, bool dbg);
    ~ZoltanMesh();
    Migration* run(MeshTag* w, double tol, int mult);
  public:
    Mesh* mesh;
    MeshTag* weights;
    bool isLocal;
    int method;
    int approach;
    bool debug;
    double tolerance;
    int multiple;
    Numbering* local;
    DynamicArray<MeshEntity*> elements;
    GlobalNumbering* global;
    MeshTag* opposite;
};

}

#endif
