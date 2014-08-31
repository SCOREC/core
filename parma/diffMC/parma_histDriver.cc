#include "parma_hist.h"
#include <assert.h>

int main() {
  //int v[] = {1,2,5,7,11};
  //int v[] = {22,1,31,3,12,102,28,27,38};
  int v[] = {22,1,31,3,12,102,28,27,38,123,110,87,65,89};
  const int vSz = sizeof(v)/sizeof(int);
  hist h;
  for(int i=0; i<vSz; i++)
    h.add(v[i]);
  assert(vSz == h.getSz());
  h.print(2);
  h.print(3);
  h.print(4);
  h.print(5);
  h.clear();
  assert(0==h.getSz());
  return 0;
}
