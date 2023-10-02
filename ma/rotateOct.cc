#include <cstdio>

int main()
{
  int starts[6][6] =
  {{0,1,2,3,4,5}
  ,{1,2,0,4,5,3}
  ,{2,3,0,1,5,4}
  ,{3,2,5,4,0,1}
  ,{4,3,5,1,0,2}
  ,{5,1,4,3,2,0}};
  for (int i=0; i < 6; ++i)
  {
    int v[6];
    for (int j=0; j < 4; ++j)
    {
      int* in = starts[i]+1;
      int* out = v+1;
      for (int k=0; k < 4; ++k)
        out[k] = in[(k+j)%4];
      v[0] = starts[i][0];
      v[5] = starts[i][5];
      if (i*4+j==0)
        printf("{");
      else
        printf(",");
      printf("{%d,%d,%d,%d,%d,%d}\n",
          v[0],v[1],v[2],v[3],v[4],v[5]);
    }
  }
  printf("};\n");
  return 0;
}
