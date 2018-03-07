#include <apfMIS.h>
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <cstdlib>
#include <apfShape.h>
#include <iostream>
#include <map>
#include <iomanip>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if ( argc != 4 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <out prefix>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  
  gmi_register_mesh();

  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);

  //We will create a field with the coloring as the values on each element
  apf::Field* coloring = apf::createField(m,"colors",apf::SCALAR,
                                          apf::getConstant(m->getDimension()));

  //Set up a MIS with primary type as elements and adjacencies as vertices.
  apf::MIS* mis = apf::initializeMIS(m,m->getDimension(),0);

  //Map to track the colors and size of the associated independent set
  std::map<int,int> colors;

  while (apf::getIndependentSet(mis)) {
    //This for loop can be thread parallized safetly
    for (int i=0;i<mis->n;i++) {
      //Independent work can be done here
      apf::setScalar(coloring,mis->ents[i],0,mis->color);
      ++colors[mis->color];
    }
  }
  apf::finalizeMIS(mis); 

  //A second coloring over the vertices
  apf::Field* coloring2 = apf::createField(m,"colors2",apf::SCALAR,
                                           m->getShape());
  //Another MIS example where primary types are vertices and
  //    adjacencies are edges
  mis = apf::initializeMIS(m,0,1);

  //second map to track colors_vertex_edge
  std::map<int,int> colors2;

  while (apf::getIndependentSet(mis)) {
    for (int i=0;i<mis->n;i++) {
      apf::setScalar(coloring2,mis->ents[i],0,mis->color);
      ++colors2[mis->color];
    }
  }
  apf::finalizeMIS(mis);

  std::cout << "MIS tests for " << argv[1] << ", " << argv[2] << std::endl;
  std::cout << "Test 1:" << std::endl;
  for (std::map<int,int>::iterator itr = colors.begin(); itr != colors.end(); ++itr) {
    std::cout << std::setw(6) << itr->second;
  }
  std::cout << std::endl;
  std::cout << "Test 2:" <<std::endl;
  for (std::map<int,int>::iterator itr = colors2.begin(); itr != colors2.end(); ++itr) {
    std::cout << std::setw(6) << itr->second;
  }
  std::cout << std::endl;

  //Vtk files should have a coloring of the elements
  //    to prove the coloring forms independent sets
  apf::writeVtkFiles(argv[3], m);
  m->destroyNative();
  apf::destroyMesh(m);

  PCU_Comm_Free();
  MPI_Finalize();
}

