#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <lionPrint.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#include <cstdlib> //exit and exit_failure
#include <PCU.h>
#include <mpi.h>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  if ( argc != 2 ) {
    if ( !PCUObj.Self() )
      printf("Usage: %s <model>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();

  gmi_model* g = 0;
  g = gmi_load(argv[1]);

  const char* tf[2] = {"false", "true"};
  printf("supports parametric evaluations: %s\n", tf[gmi_can_eval(g)]);

  for(int dim=3; dim>1; dim--) {
    struct gmi_iter* it = gmi_begin(g,dim);
    struct gmi_ent* ent;
    while( (ent = gmi_next(g,it)) ) {
      int tag = gmi_tag(g,ent);
      printf("dim %d tag %d\n", dim, tag);
      if( dim==2 && gmi_can_eval(g) ) {
        double minPt[2]; double maxPt[2];
        double r[2];
        gmi_range(g,ent,0,r);
        minPt[0] = r[0]; maxPt[0] = r[1];
        gmi_range(g,ent,1,r);
        minPt[1] = r[0]; maxPt[1] = r[1];
        double pos[3];
        gmi_eval(g,ent,minPt,pos);
        printf("  min %.3f %.3f %.3f\n", pos[0], pos[1], pos[2]);
        gmi_eval(g,ent,maxPt,pos);
        printf("  max %.3f %.3f %.3f\n", pos[0], pos[1], pos[2]);
      }
      struct gmi_set* down = gmi_adjacent(g,ent,dim-1);
      printf("  %d bounding ents with dimension %d ", down->n, dim-1);
      for(int i=0; i<down->n; i++) {
        printf(" %d ", gmi_tag(g,down->e[i]));
      }
      printf("\n");
      gmi_free_set(down);
    }
    gmi_end(g,it);
    printf("-----------\n");
  }

  gmi_destroy(g);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  }
  MPI_Finalize();
}

