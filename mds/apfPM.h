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

#ifndef MDS_PART_H
#define MDS_PART_H

#include <apfMesh2.h>

namespace apf {

struct PME
{
  PME(Parts const& i)
  {
    owner = -1;
    ids.assign(i.begin(), i.end());
    refs = 0;
  }
  bool operator<(PME const& other) const
  {
    return ids < other.ids;
  }
  int owner;
  std::vector<int> ids;
  int refs;
};

typedef std::set<PME> PM;

PME* getPME(PM& ps, Parts const& ids);
void putPME(PM& ps, PM* p);
void updateOwners(Mesh* m, PM& ps);

void initResidence(Mesh2* m, int dim);
void stitchMesh(Mesh2* m);

}

#endif
