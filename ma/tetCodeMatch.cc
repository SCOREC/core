#include "maTables.h"
#include <cstdlib>
#include <cstdio>

/* This program automatically generates the
   matching table from all edge configurations
   into the canonical configurations.
   It uses a brute force approach at many levels,
   but runs quickly once.
   This can be seen as a performance optimization
   compared to doing this search during meshadapt. */

using namespace ma;

int getOldEdge(int newEdge, int rotation)
{
  int const* new_to_old_verts = tet_rotation[rotation];
  int const* new_verts = tet_edge_verts[newEdge];
  int old_verts[2];
  old_verts[0] = new_to_old_verts[new_verts[0]];
  old_verts[1] = new_to_old_verts[new_verts[1]];
  for (int i=0; i < 6; ++i)
    if ((tet_edge_verts[i][0] == old_verts[0] &&
         tet_edge_verts[i][1] == old_verts[1])||
        (tet_edge_verts[i][0] == old_verts[1] &&
         tet_edge_verts[i][1] == old_verts[0]))
      return i;
  fprintf(stderr,"edge %d does not map under rotation %d\n",
      newEdge,rotation);
  abort();
}

int getNewCode(int oldCode, int rotation)
{
  int newCode = 0;
  for (int newEdge=0; newEdge < 6; ++newEdge)
  {
    int oldEdge = getOldEdge(newEdge,rotation);
    if (oldCode & (1<<oldEdge))
      newCode |= (1<<newEdge);
  }
  return newCode;
}

void findMatch(int code, int& rotation, int& index)
{
  for (rotation = 0; rotation < 12; ++rotation)
  {
    int newCode = getNewCode(code,rotation);
    for (int i=0; i < tet_edge_code_count; ++i)
      if (newCode == tet_edge_codes[i])
      {index=i; return;}
  }
  fprintf(stderr,"code 0x%02X has no canonical equivalent\n",code);
}

int main()
{
  printf("CodeMatch const tet_code_match[(1<<6)] =\n");
  for (int code = 0; code < (1<<6); ++code)
  {
    int rotation, index;
    findMatch(code,rotation,index);
    if (code==0) printf("{");
    else printf(",");
    printf("{%2d,%2d}\n",rotation,index);
  }
  printf("};\n");
}
