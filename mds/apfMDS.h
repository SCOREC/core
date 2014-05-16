/*
   Copyright 2014 Dan Ibanez

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#ifndef APFMDS_H
#define APFMDS_H

struct gmi_model;

namespace apf {

class Mesh;
class Mesh2;
class MeshTag;
class MeshEntity;

Mesh2* createMdsMesh(gmi_model* model, Mesh* from);

Mesh2* loadMdsMesh(const char* modelfile, const char* meshfile);

MeshTag* numberMdsMesh(Mesh2* mesh);

void defragMdsMesh(Mesh2* mesh);

int getMdsId(MeshEntity* e);

gmi_model* getMdsModel(Mesh2* mesh);

}

#endif
