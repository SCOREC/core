#include <phSnap.h>
#include <cstdio>

#ifdef __cplusplus
extern"C"{
#endif
void sim_get_pos_on_surf (double dx, double dy, double dz, int id) {
  printf("This is a dummy function. Please compile with core at phastaChef level.\n");
  printf("coord=(%f,%f,%f); ID=%d\n",dx,dy,dz,id);
}

#ifdef __cplusplus
}
#endif

