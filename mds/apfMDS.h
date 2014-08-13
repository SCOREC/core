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

int getMdsIndex(MeshEntity* e);

MeshEntity* getMdsEntity(int apfType, int index);

void splitMdsMesh(Mesh2* m, Migration* plan, int n, void (*runAfter)(Mesh2*));

bool alignMdsMatches(Mesh2* in);

void deriveMdsModel(Mesh2* in);

void changeMdsDimension(Mesh2* in, int d);

}

#endif
