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

Mesh2* makeEmptyMdsMesh(gmi_model* model, int dim, bool isMatched);

Mesh2* addMdsPart(Mesh2* original);

Mesh2* loadMdsMesh(const char* modelfile, const char* meshfile);

Mesh2* createMdsMesh(gmi_model* model, Mesh* from);

void defragMdsMesh(Mesh2* mesh);

int getMdsId(MeshEntity* e);

gmi_model* getMdsModel(Mesh2* mesh);

}

#endif
