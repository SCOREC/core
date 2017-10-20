#ifndef PH_SNAP_H
#define PH_SNAP_H

#include <apf.h>
#include <apfMesh.h>
#include<stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

void core_measure_mesh (double x1[], double x2[], double x3[], int numnp,
                        double minq);

#ifdef __cplusplus
}
#endif

#endif
