#ifndef __APF_MIS__
#define __APF_MIS__
#include "apfMesh.h"

namespace apf{  
  class MIS {
  public:
    MIS(Mesh* mesh, int vtx_dim_, int edge_dim_);
    Mesh* m;
    const int vtx_dim;
    const int edge_dim;
    MeshEntity** ents;
    int n;
    int color;
  };

  //Intializes the independent set structure
  MIS* initializeMIS(Mesh* mesh, int vtx_dim, int edge_dim);
  //Gets the next independent set of elements
  bool getIndependentSet(MIS* mis);
  //Cleans up he independent set structure
  void finalizeMIS(MIS* mis);
}

#endif
