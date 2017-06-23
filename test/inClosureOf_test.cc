#include "ma.h"
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <PCU.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <SimModel.h>
#include <MeshSim.h>
#endif
#include <pcu_util.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc==2);
  const char* modelFile = argv[1];

  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();

  gmi_model* model = gmi_load(modelFile); // read the model here


  // get one of the model faces first
  gmi_ent* gf;
  gmi_iter* gi = gmi_begin(model, 2);
  gf = gmi_next(model, gi);
  gmi_end(model, gi);

  gmi_ent* g;
  // go through the model vertexes and check if they are in closure of gf
  printf("checking verts against face with tag %d\n",
      gmi_tag(model, gf));
  gi = gmi_begin(model, 0);
  while( (g = gmi_next(model, gi)) ){
    int res = gmi_is_in_closure_of(model, g, gf);
    if (res)
      printf("vertex with tag %d IS inside the face with tag %d\n",
      	  gmi_tag(model, g), gmi_tag(model, gf));
    else
      printf("vertex with tag %d IS NOT inside the face with tag %d\n",
      	  gmi_tag(model, g), gmi_tag(model, gf));
  }
  gmi_end(model, gi); // end the iterator

  // go through the model edges and check if they are in closure of gf
  printf("checking edges against face with tag %d\n",
      gmi_tag(model, gf));
  gi = gmi_begin(model, 1);
  while( (g = gmi_next(model, gi)) ){
    int res = gmi_is_in_closure_of(model, g, gf);
    if (res)
      printf("edge with tag %d IS inside the face with tag %d\n",
      	  gmi_tag(model, g), gmi_tag(model, gf));
    else
      printf("edge with tag %d IS NOT inside the face with tag %d\n",
      	  gmi_tag(model, g), gmi_tag(model, gf));
  }
  gmi_end(model, gi); // end the iterator

  // go through the model faces and check if they are in closure of gf
  printf("checking faces against face with tag %d\n",
      gmi_tag(model, gf));
  gi = gmi_begin(model, 2);
  while( (g = gmi_next(model, gi)) ){
    int res = gmi_is_in_closure_of(model, g, gf);
    if (res)
      printf("face with tag %d IS inside the face with tag %d\n",
      	  gmi_tag(model, g), gmi_tag(model, gf));
    else
      printf("face with tag %d IS NOT inside the face with tag %d\n",
      	  gmi_tag(model, g), gmi_tag(model, gf));
  }
  gmi_end(model, gi); // end the iterator

  gmi_destroy(model); // deleting the model

  PCU_Comm_Free();
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  MPI_Finalize();
}

