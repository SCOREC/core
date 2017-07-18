#ifndef PH_SNAP_H
#define PH_SNAP_H

#include <apf.h>
#include <apfMesh.h>
#include<stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

void pass_mesh_to_phasta(apf::Mesh* mesh);
void sim_get_pos_on_surf (double dx[], double dy[], double dz[], int numnp,
                          double px[], double py[], double pz[]);

#ifdef __cplusplus
}
#endif

#endif
