#include "maTables.h"
#include <cstdio>

struct Case
{
  int code;
  ma::CodeMatch match;
};

Case cases[5] =
{
 {0x00//0000.0000 not split (tetrahedronize)
 ,{0,0}}
,{0x05//0000.0101 split two quad edges
 ,{0,1}}
,{0x0A//0000.1010 split two quad edges
 ,{1,1}}
,{0xFF//1111.1111 split all edges
 ,{0,2}}
};

int main()
{
  for (int i=0; i < (1<<8); ++i)
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
