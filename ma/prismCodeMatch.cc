#include "maTables.h"
#include <apfMesh.h>
#include <cstdlib>
#include <cstdio>

int getOldEdge(int newEdge, int rotation)
{
  int const* new_to_old_verts = ma::prism_rotation[rotation];
  int const* new_verts = apf::prism_edge_verts[newEdge];
  int old_verts[2];
  old_verts[0] = new_to_old_verts[new_verts[0]];
  old_verts[1] = new_to_old_verts[new_verts[1]];
  for (int i=0; i < 9; ++i)
    if ((apf::prism_edge_verts[i][0] == old_verts[0] &&
         apf::prism_edge_verts[i][1] == old_verts[1])||
        (apf::prism_edge_verts[i][0] == old_verts[1] &&
         apf::prism_edge_verts[i][1] == old_verts[0]))
      return i;
  fprintf(stderr,"edge %d does not map under rotation %d\n",newEdge,rotation);
  fprintf(stderr,"new edge %d\n",newEdge);
  fprintf(stderr,"new edge verts %d %d\n", new_verts[0], new_verts[1]);
  fprintf(stderr,"old edge verts %d %d\n", old_verts[0], old_verts[1]);
  abort();
}

int getNewCode(int oldCode, int rotation)
{
  int newCode = 0;
  for (int newEdge=0; newEdge < 9; ++newEdge)
  {
    int oldEdge = getOldEdge(newEdge,rotation);
    if (oldCode & (1<<oldEdge))
      newCode |= (1<<newEdge);
  }
  if (rotation == 0)
    assert(oldCode == newCode);
  return newCode;
}

void findMatch(int code, int& rotation, int& index)
{
  for (rotation = 0; rotation < 6; ++rotation)
  {
    int newCode = getNewCode(code,rotation);
    for (int i=0; i < ma::prism_edge_code_count; ++i) {
      if (newCode == ma::prism_edge_codes[i]) {
        index = i;
        return;
      }
      if (code == ma::prism_edge_codes[i])
        fprintf(stderr,"prism edge code %d not matching itself?\n",code);
    }
  }
//fprintf(stderr,"code 0x%02X has no canonical equivalent\n",code);
  rotation = 0;
  index = -1;
}

int main()
{
  printf("CodeMatch const prism_code_match[(1<<9)] =\n");
  for (int code = 0; code < (1<<9); ++code)
  {
    int rotation, index;
    findMatch(code,rotation,index);
    if (code==0) printf("{");
    else printf(",");
    printf("{%d,%d}\n",rotation,index);
  }
  printf("};\n");
  for (int i=0; i < ma::prism_edge_code_count; ++i)
    fprintf(stderr,"canonical code %d is %d\n",i,ma::prism_edge_codes[i]);
}

