#include "maTables.h"
#include <cstdlib>
#include <cstdio>

/* This program automatically generates the
   matching table from the projection code
   of a vertex onto its opposite face into
   the rotation and template used to correct the
   sliver tet.

   There are two types of sliver tets. the first
   has two nearly touching edges, the second has
   a vertex nearly touching a face. */

/*
   The codes are as follows depending on
   the projections of v3 onto triangle 0-1-2:
              |            /
       001=1  |  101=5    /  100=4
              |          /
   ---------v0+-------v2+-----------
              |        /
       011=3  | 111=7 /
              |      /
              |     /
              |    /  110=6
              |   /
              |  /
              | /
             v1/
              +
             /|
            / |
           /  |
          /   |
         /    |
        /010=2|
*/

/* the canonical orientation for an edge-edge
   sliver will be as follows:

   3+------+2
    |\    /|
    | \  / |
    |  \/  |
    |  /\  |
    | /  \ |
    |/    \|
   0+------+1

   such that edges (0,2) and (1,3) are the
   "nearly touching" edges.
*/

/* the canonical orientation for a vertex-face
   sliver will be such that vertex 3 is nearly
   touching triangle (0,1,2) */

enum
{
  edge_edge=0,
  vert_face=1
};

int const codeToType[8]=
{-1,//0
 vert_face,//1
 vert_face,//2
 edge_edge,//3
 vert_face,//4
 edge_edge,//5
 edge_edge,//6
 vert_face //7
};

/* for the vertex-face cases,
   indicates the vertex */
int const codeToVert[8]=
{-1,//0
 0, //1
 1, //2
 -1,//3
 2, //4
 -1,//5
 -1,//6
 3  //7
};

/* for the edge-edge cases, indicates the edges
   as vertex pairs */
int const codeToEdges[8][2][2]=
{{{-1,-1},{-1,-1}},//0
 {{-1,-1},{-1,-1}},//1
 {{-1,-1},{-1,-1}},//2
 {{ 0, 1},{ 2, 3}},//3
 {{-1,-1},{-1,-1}},//4
 {{ 0, 2},{ 1, 3}},//5
 {{ 1, 2},{ 0, 3}},//6
 {{-1,-1},{-1,-1}} //7
};

int findVertFaceRotation(int vert)
{
  for (int r=0; r < 12; ++r)
  {
    int const* v = ma::tet_rotation[r];
    if (v[3]==vert) return r;
  }
  assert(!"vert-face rotation not found");
}

bool isSameEdge(int const a[2], int const b[2])
{
  if ((a[0]==b[0])&&(a[1]==b[1])) return true;
  if ((a[0]==b[1])&&(a[1]==b[0])) return true;
  return false;
}

bool isSameEdgePair(int const a[2][2], int const b[2][2])
{
  if (isSameEdge(a[0],b[0])&&isSameEdge(a[1],b[1])) return true;
  if (isSameEdge(a[0],b[1])&&isSameEdge(a[1],b[0])) return true;
  return false;
}

void rotateEdge(int const a[2], int r, int b[2])
{
  int const* newToOldVerts = ma::tet_rotation[r];
  int oldToNewVerts[4];
  for (int i=0; i < 4; ++i) oldToNewVerts[newToOldVerts[i]]=i;
  b[0] = oldToNewVerts[a[0]];
  b[1] = oldToNewVerts[a[1]];
}

void rotateEdgePair(int const a[2][2], int r, int b[2][2])
{
  rotateEdge(a[0],r,b[0]);
  rotateEdge(a[1],r,b[1]);
}

int findEdgeEdgeRotation(int const edges[2][2])
{
  int const canonicalEdges[2][2]={{0,2},{1,3}};
  int rotatedEdges[2][2];
  for (int r=0; r < 12; ++r)
  {
    rotateEdgePair(edges,r,rotatedEdges);
    if (isSameEdgePair(rotatedEdges,canonicalEdges))
      return r;
  }
  assert(!"edge-edge rotation not found");
}

int main()
{
  printf("CodeMatch const sliver_code_match[8] =\n");
  printf("{{-1,-1}\n");
  for (int code = 1; code < 8; ++code)
  {
    int type = codeToType[code];
    int rotation;
    if (type==edge_edge) rotation = findEdgeEdgeRotation(codeToEdges[code]);
    if (type==vert_face) rotation = findVertFaceRotation(codeToVert[code]);
    printf(",{%d,%d}\n",rotation,type);
  }
  printf("};\n");
}

