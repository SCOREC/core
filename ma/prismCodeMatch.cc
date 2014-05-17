#include "maTables.h"
#include <cstdio>

struct Case
{
  int code;
  ma::CodeMatch match;
};

Case cases[5] =
{
 {0x000//000.000.000 not split (tetrahedronize)
 ,{0,0}}
,{0x041//001.000.001 split two horizontal edges (0)
 ,{0,1}}
,{0x082//010.000.010 split two horizontal edges (1)
 ,{1,1}}
,{0x104//100.000.100 split two horizontal edges (2)
 ,{2,1}}
,{0x1FF//111.111.111 split all edges
 ,{0,2}}
};

int main()
{
  for (int i=0; i < (1<<9); ++i)
  {
    for (int j=0; j < 5; ++j)
      if (i==cases[j].code)
      {
        printf(",{%d,%d}\n",
            cases[j].match.rotation,
            cases[j].match.code_index);
        goto nextcode;
      }
    printf(",{0,-1}\n");
nextcode:
    continue;
  }
}
