/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#ifndef APFMDS_H
#define APFMDS_H

struct gmi_model;

namespace apf {

class Mesh;
class Mesh2;
class MeshTag;
class MeshEntity;
class Migration;

Mesh2* makeEmptyMdsMesh(gmi_model* model, int dim, bool isMatched);

Mesh2* loadMdsMesh(const char* modelfile, const char* meshfile);
Mesh2* loadMdsMesh(gmi_model* model, const char* meshfile);

Mesh2* createMdsMesh(gmi_model* model, Mesh* from);

void reorderMdsMesh(Mesh2* mesh);

void splitMdsMesh(Mesh2* m, Migration* plan, int n, void (*runAfter)(Mesh2*));

bool alignMdsMatches(Mesh2* in);
bool alignMdsRemotes(Mesh2* in);

void deriveMdsModel(Mesh2* in);

void changeMdsDimension(Mesh2* in, int d);

/** \brief returns the dimension-unique index for this entity */
int getMdsIndex(Mesh2* in, MeshEntity* e);

/** \brief retrieve an entity by dimension and index
  \details indices follow iteration order, so this
  function is equivalent to iterating (index) times,
  but is actually much faster than that. */
MeshEntity* getMdsEntity(Mesh2* in, int dimension, int index);

}

#endif
