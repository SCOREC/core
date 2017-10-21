#include <phastaChef.h>
#include <gmi_mesh.h>

#include "apfSIM.h"
#include "gmi_sim.h"
#include <SimUtil.h>
#include <SimModel.h>
#include <SimPartitionedMesh.h>
#include <SimMeshTools.h>

#include <pcu_util.h>
#include <cstdio>

#ifdef __cplusplus
extern"C"{
#endif

void core_measure_mesh (double x1[], double x2[], double x3[], int numnp,
                        double minq) {
  (void) x1;
  (void) x2;
  (void) x3;
  (void) numnp;
  (void) minq;
  printf("This is a dummy function. Please compile with core at phastaChef level.\n");
}

#ifdef __cplusplus
}
#endif

